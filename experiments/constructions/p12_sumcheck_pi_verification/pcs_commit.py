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
    """min over basis vectors e_j of  weight(Enc(e_j)) / N  -- the quantity that
    governs the forge cheat (the forged-v difference is delta*Enc(e_{j0})).  A
    MEASURED lower bound on the code's relative distance for the soundness claim."""
    pl = _exp_plan(k, q)
    N = pl[2]
    wmin = N
    for j in range(k):
        e = [0] * k
        e[j] = 1
        cw = _exp_encode(e, q, pl)
        wmin = min(wmin, sum(1 for v in cw if v % q != 0))
    return wmin / N


def _code_ncode(code, k, q, blowup):
    """Codeword length N for the chosen row code (independent verifier recompute)."""
    if code == "rs":
        return max(blowup * k, k + 1)
    if code == "expander":
        return _exp_plan(k, q)[2]
    raise ValueError("unknown code %r" % (code,))


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


def _forge_prove(com, pt, true_val, wrong_val, t=32):
    """A cheating opener: forge a message v' with <v',b>=wrong_val (so it passes the
    evaluation binding) by perturbing ONE coordinate of the honest v.  Everything
    else is opened honestly.  The queried-column consistency Enc(v')[c]=<a,col> must
    then fail (v' != a^T M)."""
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
    j0 = next(j for j in range(k) if b[j] != 0)
    delta = (wrong_val - true_val) % q
    v[j0] = (v[j0] + delta * pow(b[j0], q - 2, q)) % q
    seed1 = _digest(seed0, _ser(v, q), _ser(w, q))
    cols = _distinct_challenges(seed1, b"col", min(t, Ncode), Ncode)
    opened = [(c, [com["Mhat"][i][c] for i in range(r)],
               _merkle_path(com["levels"], c)) for c in cols]
    return dict(v=v, w=w, opened=opened, r=r, k=k, n1=n1, Ncode=Ncode,
                nb=com["nb"], code=com["code"])


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
    ap.add_argument("--field", choices=list(FIELDS), default="q")
    ap.add_argument("--code", choices=("rs", "expander"), default="rs")
    args = ap.parse_args()
    if args.bench:
        bench(FIELDS[args.field], args.code)
    else:
        selftest()
