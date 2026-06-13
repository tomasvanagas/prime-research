#!/usr/bin/env python3
"""
pcs_commit.py -- S508, a hash-based multilinear polynomial commitment (tensor /
Ligero-Brakedown style) with a SUB-sqrt(x) verifier, for the chain's two
one-time S_0 base closes.

WHERE THIS SITS.  The compressed pi(x) chain (compressed_layer.run_chain) is now
per-layer LEAF-EVAL-FREE (S506 batch_ub): every per-layer carried claim threads
to the base, and the only verifier work whose cost grows with x is ONE-TIME --
the two S_0 base closes, each a DIRECT mle_eval over an O(sqrt x)=2^nb-size table
(`C == mle_eval(S0, z)` at compressed_layer.py:1047).  The whole-chain verifier
op-count (S507 `--bench-verifier-ops`) is therefore Theta(sqrt x): per-layer 0,
one-time 2*2^nb.  That one-time 2^nb term is the LAST thing between the current
state and a SUCCINCT (polylog / sub-sqrt x) verifier.

WHAT THIS BUILDS.  A polynomial commitment for a multilinear S~ over nb variables
whose EVALUATION-PROOF VERIFIER is O(sqrt(2^nb) * polylog) = O(x^{1/4} * polylog),
i.e. sub-sqrt(x).  Tensor / Ligero-Brakedown construction:

  * Reshape the 2^nb-entry table S (MSB-first) into an r x k matrix M, r=2^{n1},
    k=2^{n2}, n1+n2=nb (row bits = high n1, col bits = low n2: M = S.reshape(r,k)).
  * eq~(pt, .) factors over the bit split, eq~(pt,w)=eq~(pt_hi,w_hi)*eq~(pt_lo,w_lo),
    so with the tensor vectors a=eq_table(pt_hi) (len r), b=eq_table(pt_lo) (len k):
        S~(pt) = sum_{i,j} M[i,j] a[i] b[j] = a^T M b.
  * COMMIT: Reed-Solomon-encode each ROW of M (a degree-<k message) to a length-N
    codeword (N=blowup*k evaluation points), giving Mhat (r x N); Merkle-commit the
    N COLUMNS of Mhat (sha256, the collision-resistant-hash / random-oracle
    stand-in).  The commitment is the Merkle root.
  * OPEN (Fiat-Shamir, non-interactive): derive a proximity vector rho<-RO(root,pt,
    claimed) (len r); the prover sends the two combined MESSAGES v=a^T M (the
    evaluation combination, len k) and w=rho^T M (the proximity combination, len k),
    then opens t columns Q<-RO(root,pt,claimed,v,w) with Merkle paths.
  * VERIFY:  (1) evaluation binding  <v, b> == claimed;  (2) for each queried col c:
    the revealed column Mhat[:,c] hashes to the committed root (Merkle), AND
    Enc(v)[c] == <a, Mhat[:,c]>  AND  Enc(w)[c] == <rho, Mhat[:,c]>  (encoding is
    linear, so for an honest codeword matrix Enc(a^T M)[c] = a^T Mhat[:,c]).

  Verifier cost: reading v,w O(k); <v,b> O(k); per queried column two length-r dot
  products + two length-k encodings = O(t(r+k)) field ops + O(t r) hashed elements
  for the Merkle paths.  With r,k ~ 2^{nb/2} ~ x^{1/4} and t fixed, that is
  O(x^{1/4} * polylog) -- SUB-sqrt(x).  The full 2^nb table never touches the
  verifier; it lives only in the prover's commit/open.

SOUNDNESS (standard Ligero/Brakedown, under a collision-resistant hash / RO).
The Merkle root binds Mhat (collision resistance).  The random proximity test
(rho, t columns) forces the committed Mhat to be close to a matrix whose rows are
codewords, so it decodes to a unique M.  The evaluation columns then pin v=a^T M,
and <v,b>=claimed forces claimed = S~(pt) for that committed M.  A wrong claimed
value fails (1) outright (the honest v has <v,b>=true value); a prover who instead
forges v' with <v',b>=wrong has Enc(v') differing from a^T Mhat in >= N-(k-1)
columns (two degree-<k polys agree on < k points), caught by some queried column
with prob >= 1-((k-1)/N)^t; a tampered (non-codeword) committed row is caught by
the proximity column checks with the same bound; a tampered REVEALED column fails
its Merkle path.  This is a soundness ARGUMENT verified here empirically by the
cheating-prover tests, not a machine-checked proof.

HONEST SCOPE.  (i) The win is in the SCALING (op count slope ~0.25 in x, sub-sqrt),
NOT the absolute count at reachable n: the verifier touches ~ t*(2r+2k) ~ 4t*x^{1/4}
field elements, which only drops below the direct 2^nb close once 2^{nb/2} > ~4t,
i.e. nb > ~2*log2(4t) (~nb>=15, n>=30 at t=32) -- beyond the demo's reachable n,
exactly like the S506 wall-clock note.  The bench reports the slope to show
sub-sqrt scaling regardless of the crossover constant.  (ii) This trades
UNCONDITIONALITY for full succinctness: soundness now rests on the hash (RO/CRH),
whereas the rest of the chain verifier is unconditional and Õ(sqrt x).  (iii) RS
needs N <= q distinct points; fine at demo n, and a production system swaps in a
field-size-free linear code (Brakedown's expander code) behind the SAME
commit/open/verify interface -- the tensor reduction and the cost are unchanged.

WHAT WOULD FALSIFY THIS.  verify disagreeing with mle_eval on an honest opening;
a forged opening / wrong claimed value / tampered codeword row accepted (above the
field bound); the verifier op count not sub-sqrt-x (slope -> 0.5 not -> 0.25); or
threading it into run_chain(--commit-base) changing the chain verdict / claimed pi.

Usage:
  python3 pcs_commit.py --selftest
  python3 pcs_commit.py --bench
  python3 pcs_commit.py --bench --field big
"""

import argparse
import hashlib
import time

import numpy as np

from compressed_prover_mult_trace import (DEFAULT_Q as Q, BIG_Q, SMALL_Q,
                                          FIELDS, eq_table, mle_eval)


# ----------------------------------------------------------------------
# hashing / Merkle  (sha256 = the collision-resistant-hash / random-oracle stand-in)
# ----------------------------------------------------------------------

def _digest(*parts):
    h = hashlib.sha256()
    for p in parts:
        h.update(p)
    return h.digest()


def _wid(q):
    """Fixed byte width to serialize a field element of F_q unambiguously."""
    return (int(q).bit_length() + 7) // 8 + 1


def _ser(vec, q):
    """Deterministic byte serialization of a field vector (for hashing/FS)."""
    w = _wid(q)
    return b"".join(int(int(x) % q).to_bytes(w, "big") for x in vec)


def _challenges(seed, label, n, mod):
    """n field/index challenges in [0, mod) from RO(seed || label || ctr)."""
    out = []
    i = 0
    while len(out) < n:
        d = _digest(seed, label, i.to_bytes(4, "big"))
        out.append(int.from_bytes(d, "big") % mod)
        i += 1
    return out


def _distinct_challenges(seed, label, n, mod):
    """Up to n DISTINCT indices in [0, mod) (column queries)."""
    seen, out, i = set(), [], 0
    while len(out) < n and len(seen) < mod:
        c = int.from_bytes(_digest(seed, label, i.to_bytes(4, "big")), "big") % mod
        i += 1
        if c not in seen:
            seen.add(c)
            out.append(c)
    return out


def _merkle_build(leaves):
    """Binary Merkle tree over leaf digests (padded to a power of two with a
    fixed pad leaf).  Returns the level list (levels[0]=leaves, levels[-1]=[root])."""
    n = len(leaves)
    m = 1
    while m < n:
        m <<= 1
    pad = _digest(b"pad")
    level = list(leaves) + [pad] * (m - n)
    levels = [level]
    while len(level) > 1:
        level = [_digest(level[2 * i], level[2 * i + 1])
                 for i in range(len(level) // 2)]
        levels.append(level)
    return levels


def _merkle_path(levels, idx):
    path = []
    for lvl in levels[:-1]:
        path.append(lvl[idx ^ 1])
        idx >>= 1
    return path


def _merkle_check(root, idx, leaf, path):
    h = leaf
    for sib in path:
        h = _digest(h, sib) if (idx & 1) == 0 else _digest(sib, h)
        idx >>= 1
    return h == root


# ----------------------------------------------------------------------
# Reed-Solomon encoding (the linear code; coefficient form, points 0..N-1)
# ----------------------------------------------------------------------

def _enc_point(coeffs, x, q):
    """Evaluate sum_j coeffs[j] x^j at field point x (Horner).  This is Enc(.)[c]
    for evaluation point x = c, the verifier's per-column encoding (O(k))."""
    acc = 0
    for cj in reversed(coeffs):
        acc = (acc * x + int(cj)) % q
    return acc


def _encode_rows(M, Ncode, q):
    """Mhat[i][c] = Enc(M[i])[c], the RS codeword of every row at points 0..Ncode-1."""
    return [[_enc_point(row, c, q) for c in range(Ncode)] for row in M]


def _leaves_from_Mhat(Mhat, r, Ncode, q):
    return [_digest(b"col", c.to_bytes(4, "big"),
                    _ser([Mhat[i][c] for i in range(r)], q))
            for c in range(Ncode)]


# ----------------------------------------------------------------------
# Brakedown-style linear-time EXPANDER code (the field-size-free row code)
# ----------------------------------------------------------------------
#
# WHY.  The RS row code (above) needs N = blowup*k <= q distinct evaluation
# points -- a demo constraint that caps the committable table size and forced
# the large runs onto the 61-bit BIG_Q.  Brakedown's recursive expander code is
# a LINEAR, LINEAR-TIME-ENCODABLE code with constant relative distance over ANY
# field -- no distinct-point / N<=q requirement.  It plugs in behind the SAME
# commit/prove/verify interface: only (i) how a row is encoded (commit) and (ii)
# how the combined messages v,w are re-encoded (verify) change.  The tensor
# reduction S~(pt)=a^T M b and the homomorphic column check
#     <a, Mhat[:,c]> = sum_i a[i] Enc(M[i])[c] = Enc(a^T M)[c] = Enc(v)[c]
# hold for ANY linear Enc, so the soundness argument and the O(t(r+k))=O(x^{1/4})
# verifier are unchanged; only the code (and its relative distance delta) differ.
#
# CONSTRUCTION (Spielman / Druk-Ishai recursive code, as instantiated in
# Brakedown, Golovnev-Lee-Setty-Thaler-Wahby 2021).  Systematic, recursive:
#     Enc_n(x) = x  ||  Enc_{nm}(A x)  ||  B (Enc_{nm}(A x)),   nm = ceil(n/2)
# where A (nm x n) and B (pn x ne) are SPARSE matrices over F_q and the base case
# (n <= BASE) is a dense systematic code x || R x.  A and B are COLUMN-REGULAR
# (every input coordinate is assigned to `coldeg` output rows) so EVERY input
# influences the parity -- this is what kills the weight-1 codewords a plain
# low-density generator matrix (LDGM) would have, and is the measured difference
# between rel-distance ~0.003 (row-regular, columns missed) and ~0.45 (here).
# Encode is LINEAR-TIME: each level costs O(coldeg*n) field mults and the message
# lengths shrink geometrically (sum = O(coldeg*n) = O(k) = O(x^{1/4})).
#
# The matrices are derived deterministically from a FIXED public seed + the size
# n only (column indices) and reduced mod q (values), so the prover's commit and
# the verifier's re-encode build the IDENTICAL linear map -- a fixed public code.
#
# HONEST SCOPE.  The relative distance delta is the security parameter: the
# consistency / proximity tests catch a cheating prover with prob >= 1-(1-delta)^t.
# We MEASURE delta on the demo instances (selftest asserts min low-weight-codeword
# relative weight >= 0.3; it runs ~0.42-0.65) and pick t accordingly.  Brakedown's
# *analyzed* parameters give a proven delta (~0.07) with t~148 for 100-bit
# security; the construction shape here is theirs, the constants are demo-scale and
# distance-measured rather than from the asymptotic expander bound.

_EXP_SEED = b"pcs-expander-v1"        # fixed public seed for the code's matrices
_EXP_BASE = 4                         # base-case message-length threshold
_EXP_COLDEG = 4                       # column degree of the sparse matrices A, B
_EXP_PLAN_CACHE = {}                  # (n, q) -> plan; matrices are deterministic


def _exp_val(label, i, j, q):
    """A NONZERO field element from RO(seed||label||i||j) (matrix entry)."""
    v = int.from_bytes(_digest(_EXP_SEED, label, i.to_bytes(4, "big"),
                               j.to_bytes(4, "big")), "big") % q
    return v if v != 0 else 1


def _exp_sparse_colreg(nr, nc, coldeg, q, tag):
    """A COLUMN-REGULAR sparse matrix (nr x nc): every column is placed in
    `coldeg` distinct rows with nonzero values.  Returned as row lists of
    (col, val) -- every input coordinate thus reaches >= coldeg output rows."""
    rows = [[] for _ in range(nr)]
    d = min(coldeg, nr)
    for c in range(nc):
        rs = _distinct_challenges(_EXP_SEED, tag + b":r" + c.to_bytes(4, "big"),
                                  d, nr)
        for j, ri in enumerate(rs):
            rows[ri].append((c, _exp_val(tag + b":v", ri, c * 97 + j, q)))
    return rows


def _exp_dense(nr, nc, q, tag):
    """A dense matrix (nr x nc) as row lists of (col, val) -- the base parity."""
    return [[(c, _exp_val(tag + b":v", i, c, q)) for c in range(nc)]
            for i in range(nr)]


def _exp_plan(n, q):
    """Recursive code plan for message length n over F_q.  A plan is either
    ('base', n, N, R) or ('rec', n, N, A, sub_plan, B); N is the codeword length.
    Deterministic in (n, q); cached."""
    key = (n, q)
    if key in _EXP_PLAN_CACHE:
        return _EXP_PLAN_CACHE[key]
    if n <= _EXP_BASE:
        pn = n                                   # base parity length (rate 1/2)
        R = _exp_dense(pn, n, q, b"baseR%d" % n)
        plan = ("base", n, n + pn, R)
    else:
        nm = -(-n // 2)                          # ceil(n/2)
        A = _exp_sparse_colreg(nm, n, _EXP_COLDEG, q, b"A%d" % n)
        sub = _exp_plan(nm, q)
        ne = sub[2]
        pn = -(-n // 2)                          # parity length = ceil(n/2)
        B = _exp_sparse_colreg(pn, ne, _EXP_COLDEG, q, b"B%d" % n)
        plan = ("rec", n, n + ne + pn, A, sub, B)
    _EXP_PLAN_CACHE[key] = plan
    return plan


def _exp_matvec(rows, x, q, ops=None):
    """Sparse/dense mat-vec over F_q; rows are (col, val) lists.  Counts field
    multiplies into ops[0] if given (the linear-time encode op count)."""
    out = []
    cnt = 0
    for row in rows:
        acc = 0
        for (c, v) in row:
            acc += v * x[c]
            cnt += 1
        out.append(acc % q)
    if ops is not None:
        ops[0] += cnt
    return out


def _exp_encode(msg, q, plan=None, ops=None):
    """Encode a length-n message to its length-N codeword (systematic, linear).
    O(coldeg*n) = O(n) field mults; counted into ops[0] if given."""
    n = len(msg)
    if plan is None:
        plan = _exp_plan(n, q)
    if plan[0] == "base":
        _, _, _, R = plan
        return [int(m) % q for m in msg] + _exp_matvec(R, msg, q, ops)
    _, _, _, A, sub, B = plan
    y = _exp_matvec(A, msg, q, ops)
    ye = _exp_encode(y, q, sub, ops)
    z = _exp_matvec(B, ye, q, ops)
    return [int(m) % q for m in msg] + ye + z


def _exp_min_basis_rel_weight(k, q):
    """min over basis vectors e_j of  weight(Enc(e_j)) / N  for the expander code --
    the quantity that governs the forge cheat (the forged-v difference is
    delta*Enc(e_{j0})).  A MEASURED lower bound on the code's relative distance for
    the soundness claim.  Thin wrapper over the code-agnostic _min_basis_rel_weight."""
    return _min_basis_rel_weight("expander", k, q)


def _code_ncode(code, k, q, blowup):
    """Codeword length N for the chosen row code (independent verifier recompute)."""
    if code == "rs":
        return max(blowup * k, k + 1)
    if code == "expander":
        return _exp_plan(k, q)[2]
    raise ValueError("unknown code %r" % (code,))


def _rs_basis_codeword(j, N, q):
    """Enc(e_j) for the RS row code: the j-th basis message is the monomial x^j, so
    its codeword is [c^j mod q]_{c=0..N-1} (the verifier's evaluation points).
    Measured (modular exponentiation), not assumed."""
    return [pow(c, j, q) for c in range(N)]


def _rs_basis_weights(k, N, q):
    """weight(Enc(e_j)) for every basis j of the RS row code, via the INCREMENTAL
    power table c^j = c^{j-1}*c (O(N*k) field mults, measured -- no closed-form
    shortcut).  For a prime q > N-1 this gives N (j=0) and N-1 (j>=1: zero only at
    the point c=0), but we COMPUTE it rather than assume it."""
    weights = [0] * k
    for c in range(N):
        p = 1                                    # c^0
        for j in range(k):
            if p != 0:
                weights[j] += 1
            p = (p * c) % q
    return weights


def _min_basis_rel_weight(code, k, q, blowup=4, return_all=False):
    """min over basis messages e_j of weight(Enc(e_j))/N -- the conservative
    soundness parameter governing the forge cheat (the forged-v difference is
    exactly s*Enc(e_{j0}), caught at a queried column iff that column lies in
    supp(Enc(e_{j0}))).  Code-agnostic; field-size-free for the expander (its
    support is set by q-free indices).  Returns rel (or (rel, jmin, N, weights))."""
    N = _code_ncode(code, k, q, blowup)
    if code == "rs":
        if N > q:
            raise ValueError("RS needs N<=q distinct points (N=%d > q=%d)" % (N, q))
        weights = _rs_basis_weights(k, N, q)
    elif code == "expander":
        pl = _exp_plan(k, q)
        weights = []
        for j in range(k):
            e = [0] * k
            e[j] = 1
            cw = _exp_encode(e, q, pl)
            weights.append(sum(1 for v in cw if int(v) % q != 0))
    else:
        raise ValueError("unknown code %r" % (code,))
    wmin = min(weights)
    jmin = weights.index(wmin)
    rel = wmin / N
    return (rel, jmin, N, weights) if return_all else rel


def _basis_support(code, j0, k, q, blowup=4):
    """supp(Enc(e_{j0})) (the set of codeword columns the forge cheat at coordinate
    j0 perturbs) and the codeword length N.  The forge false-accepts iff none of the
    t queried columns lands in this set."""
    N = _code_ncode(code, k, q, blowup)
    if code == "rs":
        cw = _rs_basis_codeword(j0, N, q)
    elif code == "expander":
        e = [0] * k
        e[j0] = 1
        cw = _exp_encode(e, q, _exp_plan(k, q))
    else:
        raise ValueError("unknown code %r" % (code,))
    return frozenset(c for c in range(N) if int(cw[c]) % q != 0), N


# ----------------------------------------------------------------------
# the PCS:  commit / prove / verify
# ----------------------------------------------------------------------

def commit(S, q=Q, blowup=4, code="rs"):
    """Commit to the multilinear S~ given its 2^nb-value table S (MSB-first).
    Returns a commitment dict; the PUBLIC commitment is `com["root"]`, the rest
    is prover state used by prove().  Prover cost O(r*k*N)=O(x^{3/4}) (one-time,
    within the prover's Õ(x) budget).  `code` selects the linear row code:
    "rs" (Reed-Solomon, needs N<=q) or "expander" (Brakedown-style, field-size-
    free)."""
    S = [int(x) % q for x in np.asarray(S).ravel().tolist()]
    N0 = len(S)
    nb = N0.bit_length() - 1
    assert (1 << nb) == N0, "table size must be a power of two"
    n2 = nb // 2
    n1 = nb - n2
    r, k = 1 << n1, 1 << n2
    M = [[S[i * k + j] for j in range(k)] for i in range(r)]
    if code == "rs":
        Ncode = max(blowup * k, k + 1)      # distance N-k+1 >= 2
        assert Ncode <= q, "RS needs Ncode <= q distinct points (raise q for large nb)"
        Mhat = _encode_rows(M, Ncode, q)
    elif code == "expander":
        pl = _exp_plan(k, q)                # field-size-free; no N<=q requirement
        Ncode = pl[2]
        Mhat = [_exp_encode(row, q, pl) for row in M]
    else:
        raise ValueError("unknown code %r" % (code,))
    levels = _merkle_build(_leaves_from_Mhat(Mhat, r, Ncode, q))
    return dict(root=levels[-1][0], M=M, Mhat=Mhat, levels=levels,
                r=r, k=k, n1=n1, n2=n2, nb=nb, Ncode=Ncode, q=q, blowup=blowup,
                code=code)


def prove(com, pt, claimed, t=32):
    """Produce an evaluation proof that S~(pt) == claimed for the committed S."""
    q = com["q"]
    r, k, n1, Ncode = com["r"], com["k"], com["n1"], com["Ncode"]
    pt = [int(x) % q for x in pt]
    assert len(pt) == com["nb"], "point arity must match the committed table"
    a = [int(x) % q for x in np.asarray(eq_table(pt[:n1], q)).tolist()]
    M = com["M"]
    seed0 = _digest(com["root"], _ser(pt, q), _ser([claimed], q))
    rho = _challenges(seed0, b"rho", r, q)
    # combined messages (len k): v = a^T M (evaluation), w = rho^T M (proximity)
    v = [sum(a[i] * M[i][j] for i in range(r)) % q for j in range(k)]
    w = [sum(rho[i] * M[i][j] for i in range(r)) % q for j in range(k)]
    seed1 = _digest(seed0, _ser(v, q), _ser(w, q))
    cols = _distinct_challenges(seed1, b"col", min(t, Ncode), Ncode)
    opened = [(c, [com["Mhat"][i][c] for i in range(r)],
               _merkle_path(com["levels"], c)) for c in cols]
    return dict(v=v, w=w, opened=opened, r=r, k=k, n1=n1, Ncode=Ncode,
                nb=com["nb"], code=com["code"])


def verify(root, pt, claimed, proof, t=32, blowup=4, q=Q, stats=None):
    """Verify an evaluation proof.  Returns (ok, vops): ok is acceptance, vops is
    the count of field-element operations the verifier performs in its dominant
    loops (the SUB-sqrt-x leading term).  When stats is given, vops is added to
    stats['vcommit_ops'].  The full 2^nb table is never read."""
    code = proof.get("code", "rs")
    nb = len(pt)
    n2 = nb // 2
    n1 = nb - n2
    r, k = 1 << n1, 1 << n2
    Ncode = _code_ncode(code, k, q, blowup)
    vops = 0

    def done(ok):
        if stats is not None:
            stats["vcommit_ops"] = stats.get("vcommit_ops", 0) + vops
        return ok, vops

    # structural sanity (recomputed independently from nb + public params)
    if proof["r"] != r or proof["k"] != k or proof["Ncode"] != Ncode:
        return done(False)
    v = [int(x) % q for x in proof["v"]]
    w = [int(x) % q for x in proof["w"]]
    if len(v) != k or len(w) != k:
        return done(False)
    pt = [int(x) % q for x in pt]
    claimed = int(claimed) % q
    a = [int(x) % q for x in np.asarray(eq_table(pt[:n1], q)).tolist()]
    b = [int(x) % q for x in np.asarray(eq_table(pt[n1:], q)).tolist()]

    # (1) evaluation binding:  <v, b> == claimed
    ev = sum(v[j] * b[j] for j in range(k)) % q
    vops += k
    if ev != claimed:
        return done(False)

    # recompute the Fiat-Shamir challenges
    seed0 = _digest(root, _ser(pt, q), _ser([claimed], q))
    rho = _challenges(seed0, b"rho", r, q)
    seed1 = _digest(seed0, _ser(v, q), _ser(w, q))
    cols = _distinct_challenges(seed1, b"col", min(t, Ncode), Ncode)
    if [c for (c, _, _) in proof["opened"]] != cols:
        return done(False)

    # for a linear-time code the verifier re-encodes the full combined messages
    # ONCE (O(N)=O(k), the linear-time encode) and indexes the queried columns;
    # for RS it evaluates Enc(.)[c] per column (O(k) each).  Either way O(t(r+k)).
    if code == "expander":
        pl = _exp_plan(k, q)
        ops_box = [0]
        Vcw = _exp_encode(v, q, pl, ops_box)
        Wcw = _exp_encode(w, q, pl, ops_box)
        vops += ops_box[0]

    # (2) per-column: Merkle membership + linear-code consistency for v and w
    for (c, colvec, path) in proof["opened"]:
        if len(colvec) != r:
            return done(False)
        colvec = [int(x) % q for x in colvec]
        leaf = _digest(b"col", c.to_bytes(4, "big"), _ser(colvec, q))
        if not _merkle_check(root, c, leaf, path):
            return done(False)
        dot_a = sum(a[i] * colvec[i] for i in range(r)) % q
        dot_rho = sum(rho[i] * colvec[i] for i in range(r)) % q
        vops += 2 * r
        if code == "expander":
            enc_v_c, enc_w_c = Vcw[c], Wcw[c]
        else:
            enc_v_c, enc_w_c = _enc_point(v, c, q), _enc_point(w, c, q)
            vops += 2 * k
        if enc_v_c != dot_a or enc_w_c != dot_rho:
            return done(False)
    return done(True)


# ----------------------------------------------------------------------
# tests
# ----------------------------------------------------------------------

def _tamper_row(com, i0, rng):
    """Return a copy of com whose committed Mhat row i0 is corrupted to a random
    NON-codeword (Merkle root rebuilt to match -- a valid commitment to garbage),
    M left intact so prove() forms the honest v,w.  Models a prover that commits a
    matrix with a non-well-formed row (the proximity test must catch it)."""
    q, r, Ncode = com["q"], com["r"], com["Ncode"]
    Mhat = [list(row) for row in com["Mhat"]]
    Mhat[i0] = [(int(Mhat[i0][c]) + 1 + int(rng.integers(0, q - 1))) % q
                for c in range(Ncode)]
    levels = _merkle_build(_leaves_from_Mhat(Mhat, r, Ncode, q))
    com2 = dict(com)
    com2["Mhat"] = Mhat
    com2["levels"] = levels
    com2["root"] = levels[-1][0]
    return com2


def _forge_prove(com, pt, true_val, wrong_val, t=32, j0=None):
    """A cheating opener: forge a message v' with <v',b>=wrong_val (so it passes the
    evaluation binding) by perturbing ONE coordinate of the honest v.  Everything
    else is opened honestly.  The queried-column consistency Enc(v')[c]=<a,col> must
    then fail (v' != a^T M).  `j0` selects the bent coordinate (default: first with
    b[j0]!=0); the WORST-CASE adversary picks j0 = argmin_j weight(Enc(e_j)) among
    {j : b[j]!=0} -- the forge-rate Monte-Carlo passes that coordinate in."""
    q = com["q"]
    r, k, n1, Ncode = com["r"], com["k"], com["n1"], com["Ncode"]
    pt = [int(x) % q for x in pt]
    a = [int(x) % q for x in np.asarray(eq_table(pt[:n1], q)).tolist()]
    b = [int(x) % q for x in np.asarray(eq_table(pt[n1:], q)).tolist()]
    M = com["M"]
    seed0 = _digest(com["root"], _ser(pt, q), _ser([wrong_val % q], q))
    rho = _challenges(seed0, b"rho", r, q)
    v = [sum(a[i] * M[i][j] for i in range(r)) % q for j in range(k)]
    w = [sum(rho[i] * M[i][j] for i in range(r)) % q for j in range(k)]
    # bend v at a coordinate with b[j0] != 0 so <v',b> = wrong_val
    if j0 is None:
        j0 = next(j for j in range(k) if b[j] != 0)
    else:
        assert b[j0] != 0, "forge coordinate j0 must have b[j0] != 0 to pass binding"
    delta = (wrong_val - true_val) % q
    v[j0] = (v[j0] + delta * pow(int(b[j0]), q - 2, q)) % q
    seed1 = _digest(seed0, _ser(v, q), _ser(w, q))
    cols = _distinct_challenges(seed1, b"col", min(t, Ncode), Ncode)
    opened = [(c, [com["Mhat"][i][c] for i in range(r)],
               _merkle_path(com["levels"], c)) for c in cols]
    return dict(v=v, w=w, opened=opened, r=r, k=k, n1=n1, Ncode=Ncode,
                nb=com["nb"], code=com["code"])


# ----------------------------------------------------------------------
# S515 -- the distance/soundness curves (turning S514's single-point delta
# measurement into a measured curve in k, and the forge false-accept law in t)
# ----------------------------------------------------------------------

def _hyper_miss(N, W, t):
    """P(t DISTINCT random columns all miss a W-subset of [N]) = prod_{i<t}
    (N-W-i)/(N-i) -- the sampling-WITHOUT-replacement (hypergeometric) false-accept
    probability of the forge cheat.  <= (1-W/N)^t (the with-replacement bound)."""
    if W >= N:
        return 0.0
    p = 1.0
    for i in range(t):
        num = N - W - i
        if num <= 0:
            return 0.0
        p *= num / (N - i)
    return p


def distance_sweep(fields=("q", "big", "small"),
                   codes=("rs", "expander"),
                   ks=(8, 16, 32, 64, 128, 256, 512)):
    """S514 follow-on (a): turn the single-point distance claim into a MEASURED
    CURVE.  Tabulates delta = min_j weight(Enc(e_j))/N (the soundness parameter --
    forge false-accept <= (1-delta)^t) vs k for BOTH row codes over q & BIG_Q &
    SMALL_Q, up to k~512.  Confirms delta does NOT decay toward 0 as k grows.
    FALSIFIER: delta -> 0 with k (would break the bounded-distance soundness)."""
    import math
    print("S515 distance sweep -- min basis-codeword relative weight")
    print("  delta = min_j |Enc(e_j)|/N   (forge false-accept prob <= (1-delta)^t)")
    print("  FALSIFIER: delta -> 0 as k grows.\n")
    summary = {}
    for code in codes:
        print(f"--- code = {code} ---")
        print(f"{'k':>5} " + "".join(f"{('N['+f+']'):>10}" for f in fields)
              + "  | " + "".join(f"{('d['+f+']'):>9}" for f in fields))
        for k in ks:
            ns, ds = [], []
            for f in fields:
                q = FIELDS[f]
                try:
                    rel, jmin, N, _ = _min_basis_rel_weight(code, k, q, return_all=True)
                    ns.append(str(N))
                    ds.append(f"{rel:.3f}")
                    summary[(code, f, k)] = rel
                except ValueError:
                    ns.append("skip(N>q)")
                    ds.append("   -")
            print(f"{k:>5} " + "".join(f"{n:>10}" for n in ns) + "  | "
                  + "".join(f"{d:>9}" for d in ds))
        # report the trend over the largest few k for the first available field:
        # per-doubling decrements + a geometric-tail extrapolation of the floor.
        for f in fields:
            xs = [k for k in ks if (code, f, k) in summary]
            if len(xs) >= 3:
                dv = [summary[(code, f, k)] for k in xs]
                decr = [dv[i] - dv[i + 1] for i in range(len(dv) - 1)]
                print(f"    [{f}] delta: k={xs[0]} -> {dv[0]:.3f}   "
                      f"k={xs[-1]} -> {dv[-1]:.3f}   (min over k = {min(dv):.3f})")
                print(f"    [{f}] per-doubling decrements (delta_k - delta_2k): "
                      + " ".join(f"{d:+.3f}" for d in decr))
                # if the decline is decelerating geometrically, estimate the floor
                pos = [d for d in decr if d > 1e-4]
                if len(pos) >= 2 and pos[-1] < pos[0]:
                    ratios = [pos[i + 1] / pos[i] for i in range(len(pos) - 1)
                              if pos[i] > 0]
                    rho = sum(ratios) / len(ratios) if ratios else 0.0
                    if 0 < rho < 1:
                        tail = decr[-1] * rho / (1 - rho)   # sum of remaining drops
                        print(f"    [{f}] decrements shrink (mean ratio {rho:.2f}); "
                              f"geometric-tail floor estimate ~ {dv[-1] - tail:.3f} "
                              f"(MEASURED decline decelerates, not -> 0)")
                else:
                    print(f"    [{f}] non-decreasing -> bounded below by {dv[0]:.3f}")
                # practical consequence: query count t for 100-bit soundness
                # ((1-delta)^t <= 2^-100  =>  t >= 100 / -log2(1-delta)).
                def _t100(d):
                    return int(math.ceil(100.0 / -math.log2(1 - d)))
                print(f"    [{f}] t for 2^-100 forge soundness: k={xs[0]} -> "
                      f"t={_t100(dv[0])},  k={xs[-1]} -> t={_t100(dv[-1])} "
                      f"(grows only modestly as delta declines)")
                break
        print()
    return summary


def forge_rate(code="expander", field="q", nb=8, trials=8000, tmax=12, seed=2024,
               verbose=True):
    """S514 follow-on (b): the empirical forge false-accept rate vs the column count
    t, against the predicted (1-delta)^t and the exact hypergeometric.  The
    WORST-CASE adversary bends the honest v at the min-weight basis coordinate
    j0=argmin_j weight(Enc(e_j)) (subject to b[j0]!=0), so the difference codeword
    has support delta*N -- the conservative bound.  Across `trials` random openings
    (independent FS challenges via random pt + wrong value), we measure the fraction
    accepted at each t.  FALSIFIER: empirical accept rate >> (1-delta)^t would mean
    the column queries are correlated (FS not behaving as a random oracle here)."""
    q = FIELDS[field]
    rng = np.random.default_rng(seed)
    S = np.array([int(rng.integers(0, q)) for _ in range(1 << nb)], dtype=object)
    com = commit(S, q, code=code)
    r, k, n1, Ncode = com["r"], com["k"], com["n1"], com["Ncode"]
    M, root = com["M"], com["root"]
    rel, jmin, N, weights = _min_basis_rel_weight(code, k, q, return_all=True)
    supp_glob = _basis_support(code, jmin, k, q)[0]
    order = sorted(range(k), key=lambda j: weights[j])     # min-weight coords first
    supp_cache = {jmin: supp_glob}
    delta = rel
    tmax = min(tmax, Ncode)
    acc = [0] * (tmax + 1)
    for _ in range(trials):
        pt = [int(rng.integers(0, q)) for _ in range(nb)]
        a = [int(x) % q for x in np.asarray(eq_table(pt[:n1], q)).tolist()]
        b = [int(x) % q for x in np.asarray(eq_table(pt[n1:], q)).tolist()]
        # worst-case admissible coordinate (b[j0]!=0); almost always = global argmin
        j0 = jmin if b[jmin] != 0 else next(j for j in order if b[j] != 0)
        if j0 not in supp_cache:
            supp_cache[j0] = _basis_support(code, j0, k, q)[0]
        supp = supp_cache[j0]
        v = [sum(a[i] * M[i][j] for i in range(r)) % q for j in range(k)]
        true = sum(v[j] * b[j] for j in range(k)) % q
        wrong = (true + 1 + int(rng.integers(0, q - 1))) % q   # any value != true
        seed0 = _digest(root, _ser(pt, q), _ser([wrong], q))
        rho = _challenges(seed0, b"rho", r, q)
        w = [sum(rho[i] * M[i][j] for i in range(r)) % q for j in range(k)]
        vf = list(v)
        vf[j0] = (vf[j0] + (wrong - true) * pow(int(b[j0]), q - 2, q)) % q
        seed1 = _digest(seed0, _ser(vf, q), _ser(w, q))
        cols = _distinct_challenges(seed1, b"col", tmax, Ncode)
        hit = tmax + 1                          # first (1-indexed) column in supp
        for i, c in enumerate(cols):
            if c in supp:
                hit = i + 1
                break
        for t in range(1, tmax + 1):
            if hit > t:                         # none of the first t columns hit
                acc[t] += 1
    import math
    emp = [acc[t] / trials for t in range(tmax + 1)]
    if verbose:
        print(f"S515 forge-rate Monte-Carlo  code={code} field={field} nb={nb} "
              f"(r={r}, k={k}, N={Ncode}); trials={trials}")
        print(f"  worst-case j0=argmin weight: delta = |Enc(e_j0)|/N = "
              f"{len(supp_glob)}/{N} = {delta:.4f}")
        print(f"  {'t':>3} {'emp_accept':>11} {'(1-d)^t':>10} {'hypergeom':>10} "
              f"{'emp/pred':>9}")
    # log-linear fit of the empirical accept rate -> effective delta
    fit_t, fit_ly = [], []
    for t in range(1, tmax + 1):
        pred_wr = (1 - delta) ** t
        pred_hg = _hyper_miss(N, len(supp_glob), t)
        ratio = emp[t] / pred_wr if pred_wr > 0 else float("nan")
        if verbose:
            print(f"  {t:>3} {emp[t]:>11.5f} {pred_wr:>10.5f} {pred_hg:>10.5f} "
                  f"{ratio:>8.2f}x")
        if emp[t] > 0:
            fit_t.append(t)
            fit_ly.append(math.log(emp[t]))
    delta_eff = float("nan")
    if len(fit_t) >= 2:
        sl = _slope(fit_t, fit_ly)              # ln(accept) ~ sl * t
        delta_eff = 1 - math.exp(sl)
        if verbose:
            print(f"  fitted delta_eff (from d ln(accept)/dt) = {delta_eff:.4f}   "
                  f"measured delta = {delta:.4f}")
    if verbose:
        print("  (emp/pred ~ <=1 confirms uncorrelated queries; >>1 would falsify "
              "FS randomness)")
    return dict(emp=emp, delta=delta, delta_eff=delta_eff, N=N,
                W=len(supp_glob), tmax=tmax, trials=trials)


def selftest():
    rng = np.random.default_rng(0)
    # 1-5. correctness + the 4 cheat classes, over q & BIG_Q & SMALL_Q, for BOTH
    #      row codes ("rs": Reed-Solomon; "expander": Brakedown-style field-size-free).
    for code in ("rs", "expander"):
        for q in (Q, BIG_Q, SMALL_Q):
            for nb in (1, 2, 3, 4, 6, 8):
                for _ in range(5):
                    S = np.array([int(rng.integers(0, q)) for _ in range(1 << nb)],
                                 dtype=object)
                    pt = [int(rng.integers(0, q)) for _ in range(nb)]
                    true = int(mle_eval(S, pt, q)) % q
                    com = commit(S, q, code=code)
                    pf = prove(com, pt, true)

                    # 1. honest opening verifies and AGREES with mle_eval
                    ok, vops = verify(com["root"], pt, true, pf, q=q)
                    assert ok, (code, q, nb, "honest rejected")
                    assert vops > 0

                    # 2. wrong claimed value rejected (honest prover, wrong target)
                    wrong = (true + 1) % q
                    pf_w = prove(com, pt, wrong)
                    assert not verify(com["root"], pt, wrong, pf_w, q=q)[0], \
                        (code, q, nb, "wrong claim accepted")
                    # and the honest proof checked against a wrong target rejected
                    assert not verify(com["root"], pt, wrong, pf, q=q)[0], (code, q, nb)

                    # 3. forged opening (v' bent to match a wrong value) rejected
                    if nb >= 2:
                        pf_f = _forge_prove(com, pt, true, wrong)
                        # passes the evaluation binding <v',b>=wrong by construction...
                        assert (sum(int(pf_f["v"][j]) *
                                    int(np.asarray(eq_table(pt[com["n1"]:], q))
                                        .tolist()[j]) for j in range(com["k"])) % q
                                == wrong)
                        # ...but the column-consistency checks catch it
                        assert not verify(com["root"], pt, wrong, pf_f, q=q)[0], \
                            (code, q, nb, "forged opening accepted")

                    # 4. tampered (non-codeword) committed row rejected
                    if nb >= 1:
                        ct = _tamper_row(com, int(rng.integers(0, com["r"])), rng)
                        pf_t = prove(ct, pt, true)
                        assert not verify(ct["root"], pt, true, pf_t, q=q)[0], \
                            (code, q, nb, "tampered codeword row accepted")

                    # 5. tampered REVEALED column (single entry) fails its Merkle path
                    if pf["opened"]:
                        pf_c = dict(pf)
                        op = [list(o) for o in pf["opened"]]
                        op[0][1] = list(op[0][1])
                        op[0][1][0] = (int(op[0][1][0]) + 1) % q
                        pf_c["opened"] = [tuple(o) for o in op]
                        assert not verify(com["root"], pt, true, pf_c, q=q)[0], \
                            (code, q, nb, "tampered revealed column accepted")

    # 6. stats accounting: vcommit_ops accumulates the returned vops (both codes)
    q = Q
    for code in ("rs", "expander"):
        S = np.array([int(rng.integers(0, q)) for _ in range(1 << 8)], dtype=object)
        pt = [int(rng.integers(0, q)) for _ in range(8)]
        true = int(mle_eval(S, pt, q)) % q
        com = commit(S, q, code=code)
        pf = prove(com, pt, true)
        st = {}
        ok, vops = verify(com["root"], pt, true, pf, q=q, stats=st)
        assert ok and st["vcommit_ops"] == vops, (code,)

    # 7. it really is the SAME map as the chain's pad+mle_eval base close (both codes)
    import compressed_layer as _cl
    for code in ("rs", "expander"):
        for n in (8, 10, 12):
            x = (1 << n) - 1
            nb = max(1, _cl.isqrt(x).bit_length())
            _, sm, lg = _cl.compressed_lucy(x)
            for tab in (sm[0], lg[0]):
                S = _cl.pad(tab, nb)
                pt = [int(rng.integers(0, Q)) for _ in range(nb)]
                claimed = int(mle_eval(S, pt)) % Q
                com = commit(S, Q, code=code)
                ok, _ = verify(com["root"], pt, claimed, prove(com, pt, claimed), q=Q)
                assert ok, (code, n)

    # 8. MEASURED distance of the expander code: the min basis-codeword relative
    #    weight (which governs the forge cheat) is comfortably bounded away from 0
    #    -- the field-independent (column indices) structure gives the same support
    #    over q and BIG_Q, so distance is a property of the code not the field.
    for q in (Q, BIG_Q):
        for k in (8, 16, 32, 64):
            d = _exp_min_basis_rel_weight(k, q)
            assert d >= 0.3, ("expander distance too low", q, k, d)
    # field-independence of the codeword SUPPORT (the matrices' indices are q-free)
    for k in (8, 32):
        pl_q, pl_b = _exp_plan(k, Q), _exp_plan(k, BIG_Q)
        for j in range(k):
            e = [0] * k
            e[j] = 1
            supp_q = [i for i, v in enumerate(_exp_encode(e, Q, pl_q)) if v % Q]
            supp_b = [i for i, v in enumerate(_exp_encode(e, BIG_Q, pl_b)) if v % BIG_Q]
            assert supp_q == supp_b, ("expander support field-dependent", k, j)

    # 9. THE FIELD-SIZE-FREE WIN: over a TINY prime where RS cannot commit
    #    (Ncode = max(blowup*k, k+1) > q), the expander code commits, verifies,
    #    and rejects the cheats -- removing the N<=q demo constraint that capped
    #    the committable table size and forced the large runs onto a 61-bit field.
    qt = 17
    rng = np.random.default_rng(11)
    nb = 8                                          # k=16 -> RS Ncode=64 > 17
    S = np.array([int(rng.integers(0, qt)) for _ in range(1 << nb)], dtype=object)
    pt = [int(rng.integers(0, qt)) for _ in range(nb)]
    true = int(mle_eval(S, pt, qt)) % qt
    try:
        commit(S, qt, code="rs")
        raise AssertionError("RS should refuse Ncode > q")
    except AssertionError as e:
        assert "Ncode <= q" in str(e), e
    com = commit(S, qt, code="expander")
    assert com["Ncode"] > qt                        # codeword longer than the field
    assert verify(com["root"], pt, true, prove(com, pt, true), q=qt)[0]
    wrong = (true + 1) % qt
    assert not verify(com["root"], pt, wrong, prove(com, pt, wrong), q=qt)[0]
    assert not verify(com["root"], pt, wrong,
                      _forge_prove(com, pt, true, wrong), q=qt)[0]

    # 10. S515 DISTANCE SWEEP (falsifier: delta -> 0 as k grows).  The min
    #     basis-codeword relative weight is the soundness parameter; confirm it
    #     stays BOUNDED BELOW as k climbs (here to 128) for BOTH row codes over
    #     q & BIG_Q & SMALL_Q, and is field-independent for the expander (q-free
    #     indices).  RS sits near (N-1)/N; the expander floors well above 0.
    import math
    for code in ("rs", "expander"):
        for f in ("q", "big", "small"):
            q = FIELDS[f]
            ks = (8, 16, 32, 64, 128)
            ds = []
            for k in ks:
                if code == "rs" and _code_ncode("rs", k, q, 4) > q:
                    continue                            # RS: skip when N>q
                ds.append(_min_basis_rel_weight(code, k, q))
            assert len(ds) >= 3, (code, f, "too few measurable k")
            floor = 0.5 if code == "rs" else 0.40        # measured thresholds (k<=128)
            # (a) DELTA DOES NOT DECAY TO 0 -- the real falsifier.
            assert min(ds) >= floor, ("distance below floor", code, f, ds)
            # (b) the decline is DECELERATING (not collapsing): the per-doubling
            #     decrement at the largest k is no bigger than at the smallest k.
            #     RS increases (decrements <= 0); the expander declines but with a
            #     SHRINKING decrement (consistent with a positive asymptotic floor,
            #     Brakedown's proven constant relative distance).
            decr_first = ds[0] - ds[1]
            decr_last = ds[-2] - ds[-1]
            assert decr_last <= decr_first + 0.02, \
                ("delta decline accelerating (would point to ->0)", code, f, ds)
    # expander distance is field-independent (support set by q-free indices)
    for k in (16, 64):
        dq = _min_basis_rel_weight("expander", k, Q)
        db = _min_basis_rel_weight("expander", k, BIG_Q)
        dsm = _min_basis_rel_weight("expander", k, SMALL_Q)
        assert dq == db == dsm, ("expander distance field-dependent", k, dq, db, dsm)
    # the hypergeometric helper is a valid <= (1-delta)^t bound, decreasing in t
    for (N, W) in ((44, 20), (64, 63), (380, 170)):
        prev = 1.0
        for t in range(1, 9):
            hg = _hyper_miss(N, W, t)
            assert 0.0 <= hg <= (1 - W / N) ** t + 1e-12, ("hyper > with-repl", N, W, t)
            assert hg <= prev + 1e-12                    # monotone non-increasing
            prev = hg

    # 11. S515 FORGE-RATE LAW (falsifier: empirical accept >> (1-delta)^t -> the
    #     column queries are correlated).  Run the worst-case-coordinate Monte-Carlo
    #     at small nb for both codes; assert (a) the empirical accept rate is at most
    #     the with-replacement bound (1-delta)^t within binomial tolerance at every
    #     t, (b) it is monotone non-increasing in t, (c) the fitted delta_eff tracks
    #     the measured delta, (d) t=1 accept ~ (1-delta).  This is the empirical
    #     confirmation that t queries drive the false-accept to (1-delta)^t.
    for code in ("rs", "expander"):
        res = forge_rate(code=code, field="q", nb=6, trials=4000, tmax=7, seed=99,
                         verbose=False)
        emp, delta, N, W = res["emp"], res["delta"], res["N"], res["W"]
        for t in range(1, res["tmax"] + 1):
            pred = (1 - delta) ** t
            tol = 4 * math.sqrt(max(pred, 1e-9) * (1 - pred) / res["trials"]) + 0.02
            assert emp[t] <= pred + tol, \
                ("forge accept >> (1-delta)^t (correlated queries?)", code, t,
                 emp[t], pred)
            if t >= 2:
                assert emp[t] <= emp[t - 1] + 1e-9, ("accept not monotone", code, t)
        # t=1 accept ~ (1-delta) (one column catches with prob ~ delta)
        assert abs(emp[1] - (1 - delta)) <= 0.05, \
            ("t=1 accept != 1-delta", code, emp[1], 1 - delta)
        # fitted decay rate tracks the measured distance (not slower => not falsified)
        if not math.isnan(res["delta_eff"]):
            assert res["delta_eff"] >= delta - 0.10, \
                ("decay slower than (1-delta)^t", code, res["delta_eff"], delta)

    print("selftest OK")


def bench(q=Q, code="rs"):
    """The headline: the evaluation-proof VERIFIER op count is SUB-sqrt(x).  We
    tabulate, vs nb, the verifier field-op count `vcommit_ops` against the direct
    mle_eval close it replaces (O(2^nb)=O(sqrt x)), and fit the leading exponent
    in n (x=2^n-1, table size 2^nb ~ sqrt x ~ x^{1/2}, so direct ~ x^{0.5}, the
    commit verifier ~ x^{0.25}).  Falsifier: commit-verify op slope -> 0.5.  The
    expander code must keep the SAME sub-sqrt slope (it just removes the N<=q
    requirement)."""
    import math
    rng = np.random.default_rng(7)
    print(f"q = {q};  code = {code};  multilinear-eval VERIFIER op count: tensor "
          f"commitment (sub-sqrt x) vs direct mle_eval close (O(2^nb))")
    print(f"{'nb':>3} {'2^nb':>9} {'r':>5} {'k':>5} {'N':>5} {'commit_vops':>12} "
          f"{'direct_2^nb':>12} {'ratio':>7} {'proof_cols':>10}")
    rows = []
    for nb in (4, 6, 8, 10, 12, 14):
        S = np.array([int(rng.integers(0, q)) for _ in range(1 << nb)], dtype=object)
        pt = [int(rng.integers(0, q)) for _ in range(nb)]
        true = int(mle_eval(S, pt, q)) % q
        com = commit(S, q, code=code)
        pf = prove(com, pt, true)
        ok, vops = verify(com["root"], pt, true, pf, q=q)
        assert ok
        direct = 1 << nb
        rows.append((nb, vops, direct))
        print(f"{nb:>3} {1 << nb:>9} {com['r']:>5} {com['k']:>5} {com['Ncode']:>5} "
              f"{vops:>12} {direct:>12} {direct / vops:>6.2f}x "
              f"{len(pf['opened']):>10}")
    pts = rows[-4:]
    ns = [p[0] for p in pts]                      # nb increments by 2 here
    sl_commit = _slope(ns, [math.log2(p[1]) for p in pts])
    sl_direct = _slope(ns, [math.log2(p[2]) for p in pts])
    print(f"\nfitted slope d log2(ops)/d nb  (last 4 pts):  commit = {sl_commit:.3f}"
          f"  direct = {sl_direct:.3f}")
    print("(direct slope -> 1.0 per nb = Theta(2^nb)=Theta(sqrt x); commit slope "
          "-> 0.5 per nb = Theta(2^{nb/2})=Theta(x^{1/4}), i.e. SUB-sqrt x.)")
    # in chain terms nb ~ n/2, so per-n exponent alpha = slope/2: direct ~0.5, commit ~0.25
    print(f"in x: alpha_commit ~ {sl_commit / 2:.3f} (target ~0.25), "
          f"alpha_direct ~ {sl_direct / 2:.3f} (~0.5).")
    assert sl_commit < 0.75, ("commit slope not sub-sqrt", sl_commit)
    assert sl_direct > 0.9, ("direct slope not ~linear", sl_direct)


def _slope(xs, ys):
    k = len(xs)
    mx, my = sum(xs) / k, sum(ys) / k
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    den = sum((a - mx) ** 2 for a in xs)
    return num / den if den else 0.0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--bench", action="store_true")
    ap.add_argument("--distance-sweep", action="store_true",
                    help="S515: min basis-codeword rel weight delta vs k (curve)")
    ap.add_argument("--forge-rate", action="store_true",
                    help="S515: empirical forge accept rate vs t vs (1-delta)^t")
    ap.add_argument("--field", choices=list(FIELDS), default="q")
    ap.add_argument("--code", choices=("rs", "expander"), default="rs")
    ap.add_argument("--nb", type=int, default=8, help="table bits for --forge-rate")
    ap.add_argument("--trials", type=int, default=8000, help="--forge-rate trials")
    ap.add_argument("--kmax", type=int, default=512, help="largest k for --distance-sweep")
    args = ap.parse_args()
    if args.bench:
        bench(FIELDS[args.field], args.code)
    elif args.distance_sweep:
        ks = tuple(k for k in (8, 16, 32, 64, 128, 256, 512) if k <= args.kmax)
        distance_sweep(ks=ks)
    elif args.forge_rate:
        forge_rate(code=args.code, field=args.field, nb=args.nb, trials=args.trials)
    else:
        selftest()
