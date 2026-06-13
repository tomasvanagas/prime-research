#!/usr/bin/env python3
"""cert_incompressibility.py -- an information-theoretic floor probe for the
Õ(√x) certificate of π(x)=c (open item 3, lower-bound face).

THE QUESTION
------------
S509/S510 proved the chain certificate's SIZE is Θ(√x) and that this is
"layering-inherent": it is the K=π(√x) SEQUENTIAL per-prime layer reductions,
the ⌊v/p⌋-semigroup wall (open item 2). That is a statement about THESE
constructions. The complementary lower-bound question is: is a sub-√x
certificate information-theoretically IMPOSSIBLE, or merely unachieved?

A certificate's SIZE is bounded below by the INFORMATION it must convey. So we
separate two notions that are easy to conflate:

  * cert SIZE      -- one transcript block per layer × K=π(√x) layers ⇒ Θ(√x),
                      regardless of information (the construction has K layers).
  * cert INFORMATION-- the joint "hard" bits the transcript must carry that a
                      polylog verifier cannot itself predict.

If the INFORMATION is Θ(√x), NO certificate (of any kind that reconstructs the
sieve state) can be smaller than √x -- an information barrier. If it is polylog,
a sub-√x certificate is NOT ruled out by information; the √x is then purely the
construction's sequential structure, and the lower bound (if any) must be
COMPUTATIONAL (#P-hardness / the algebraic semigroup of open item 2), not
information-theoretic.

WHAT WE MEASURE
---------------
The chain pins, per layer a=1..K, the large-side count after sieving by the
first a primes -- exactly the Legendre partial-sieve function
    phi(x,a) = #{ 1<=n<=x : n has no prime factor <= p_a }.
phi(x,a) is the principled per-layer "checkpoint scalar" (the running survivor
count the chain carries; honest proxy -- see falsifiability). Against a POLYLOG
predictor pred(x,a) (the best a verifier could compute from x alone, no √x of
primes), the residual r(x,a)=phi(x,a)-pred(x,a) at INTEGER precision is the part
the certificate must actually convey. We measure its total description length

    H_chain(x) = sum_a ceil(log2(1+|r(x,a)|))        (bits to transmit residuals)

and its exponent vs x (slope in log H_chain / log x):
    ~0.5  => H_chain ~ √x  => the chain-class certificate IS information-forced.
    ~0    => H_chain polylog => √x is construction-structural; sub-√x not ruled
             out on information grounds (lower bound must be computational).

THREE predictors, increasingly generous (more is removed => smaller residual =>
a STRICTER test of info-forcing):
  (M) exact Mertens product mu=x*prod_{p<=p_a}(1-1/p)  -- PRIME-AWARE, uses √x of
      primes, so NOT a fair polylog predictor; an upper bound on what a
      prime-knowing verifier removes. Labeled accordingly.
  (A) asymptotic mu=x*e^{-gamma}/ln(p_a)               -- polylog (needs x,p_a).
  (F) degree-d polynomial fit of phi(x,a)/x vs the Buchstab variable
      u=ln x/ln p_a, d+1 = polylog params                -- the MOST generous
      polylog model; d+1<<K, so the params are negligible and the residual-of-fit
      bits are a fair lower bound on the hard information.

Plus controls:
  * single-value control: bits of phi(x,K)-related residual = π(x)-Li(x) (and
    -R(x)) ~ O(log x) -- reproduces the S36 reframing.
  * across-x π fluctuation magnitude exponent ~0.5 (the known zeta √x), to keep
    the COMPUTATION barrier (informational, single value O(log x) bits) cleanly
    distinct from the CERTIFICATION question (joint K-checkpoint info).
  * gzip of the integer residual sequence as a model-free compressibility
    cross-check on H_chain.

FALSIFIABILITY
--------------
* If H_chain(x) (predictor F, the most generous) fits exponent ~0.5 over the
  dyadic window, the chain-class √x certificate is information-forced -- a
  barrier; the construction face is then tight on information grounds too.
* If it fits ~0 (polylog), the √x is NOT information-forced; a sub-√x
  certificate is not ruled out by information, and item-3's lower bound must come
  from computational hardness (#P-hardness), not incompressibility. Either is a
  clean, falsifiable redirection of the lower-bound program.
* PROXY caveat (honest scope): phi(x,a) is the per-layer SCALAR state; the actual
  transcript also carries MLE-at-random-point hashes whose information is bounded
  ABOVE by the state's, so H_chain on phi is a faithful proxy and, if anything,
  an under-count of transcript redundancy -- it bounds the chain-CLASS certificate
  (a verifier that follows the sieve recursion), NOT every conceivable certificate
  (a universal bound needs the #P-hardness route). Stated in the results file.
* Integer precision is load-bearing (CLOSED_PATHS row 737, S65): relative
  smoothness does NOT help -- the certificate needs ±0 exact values, so residuals
  are rounded to integers before counting bits.

Run:
  python3 cert_incompressibility.py --selftest
  python3 cert_incompressibility.py --xmax 16777216
  python3 cert_incompressibility.py --xmax 16777216 --degree 8
"""
import argparse
import gzip
import math
import sys
import numpy as np

GAMMA = 0.5772156649015328606


# ----------------------------------------------------------------------------
# Exact number theory
# ----------------------------------------------------------------------------
def sieve_spf(n):
    """Smallest-prime-factor sieve up to n (inclusive). spf[k] = least prime
    dividing k for k>=2; spf[0]=spf[1]=0. Also returns the prime list."""
    spf = np.zeros(n + 1, dtype=np.int64)
    for i in range(2, n + 1):
        if spf[i] == 0:  # i is prime
            spf[i:n + 1:i][spf[i:n + 1:i] == 0] = i
    primes = np.nonzero(spf == np.arange(n + 1))[0]
    primes = primes[primes >= 2]
    return spf, primes


def phi_checkpoints(x, spf, primes):
    """phi(x,a) for a=0..K, K=π(√x). phi(x,a)=#{1<=n<=x : spf(n)>p_a} (n=1 counts).

    For n<=x, spf(n)<=√x OR n is prime (composite => spf<=√n<=√x). So a number is
    "killed at layer a" iff spf(n)=p_a with p_a<=√x; primes>√x and n=1 survive all
    K layers. Returns (phi array length K+1, K, sqrtx, primes_le_sqrtx)."""
    sqrtx = int(math.isqrt(x))
    K = int(np.searchsorted(primes, sqrtx, side='right'))  # π(√x)
    pls = primes[:K]
    spf_vals = spf[2:x + 1]
    # layer index (0-based) of spf among ALL primes; primes>√x land at idx>=K.
    layer = np.searchsorted(primes, spf_vals)  # exact: spf is itself a prime
    kill = np.bincount(layer, minlength=len(primes))[:K]  # kill[j]=#{n: spf=p_{j+1}}
    phi = np.empty(K + 1, dtype=np.int64)
    phi[0] = x
    phi[1:] = x - np.cumsum(kill)
    return phi, K, sqrtx, pls


def mertens_product(primes_le):
    """Cumulative prod_{i<=a}(1-1/p_i) for a=1..len. PRIME-AWARE (uses √x primes)."""
    return np.cumprod(1.0 - 1.0 / primes_le.astype(np.float64))


def li(x):
    """Logarithmic integral li(x)=PV∫_0^x dt/ln t via the convergent series
    li(x)=γ+ln(ln x)+Σ_{k>=1} (ln x)^k/(k·k!)."""
    L = math.log(x)
    s = 0.0
    term = 1.0
    for k in range(1, 200):
        term *= L / k          # term = L^k / k!
        add = term / k
        s += add
        if k > L and abs(add) < 1e-12 * abs(s):
            break
    return GAMMA + math.log(L) + s


def riemann_R(x, mu_small):
    """Riemann's R(x)=Σ_{n>=1} μ(n)/n · li(x^{1/n}). mu_small: Möbius array."""
    total = 0.0
    for n in range(1, len(mu_small)):
        m = mu_small[n]
        if m == 0:
            continue
        root = x ** (1.0 / n)
        if root < 2.0:
            break
        total += (m / n) * li(root)
    return total


def mobius(n):
    spf, _ = sieve_spf(n)
    mu = np.ones(n + 1, dtype=np.int64)
    mu[0] = 0
    for k in range(2, n + 1):
        m = k
        cur = mu[k]
        sign = 1
        last = 0
        ok = True
        while m > 1:
            p = int(spf[m])
            if p == last:
                ok = False
                break
            last = p
            m //= p
            sign = -sign
        mu[k] = sign if ok else 0
    return mu


# ----------------------------------------------------------------------------
# Information / compressibility metrics
# ----------------------------------------------------------------------------
def magbits(resid):
    """Σ ceil(log2(1+|r|)) -- bits to transmit each integer residual if they were
    INDEPENDENT (no cross-layer redundancy). This is an UPPER bound and, because
    there are K=π(√x) layers each typically contributing >=1 bit, is ~√x almost by
    construction -- it is the reference, NOT the discriminator."""
    r = np.rint(resid).astype(np.int64)
    a = np.abs(r).astype(np.float64)
    return float(np.sum(np.where(a < 0.5, 0.0, np.ceil(np.log2(1.0 + a)))))


def ar_residual_bits(resid, order=3):
    """Conditional bits: fit a linear AR(order) predictor of r(a) from its
    `order` predecessors and report magbits of the AR residual. If the K layers
    are jointly redundant (each predictable from neighbours), AR bits << magbits;
    if they are independent surprises, AR bits ~ magbits."""
    r = np.rint(resid).astype(np.float64)
    n = len(r)
    if n <= order + 2:
        return magbits(resid)
    X = np.column_stack([r[order - 1 - j: n - 1 - j] for j in range(order)])
    y = r[order:]
    coef, *_ = np.linalg.lstsq(np.column_stack([X, np.ones(len(y))]), y, rcond=None)
    pred = np.column_stack([X, np.ones(len(y))]) @ coef
    return magbits(y - pred)


def gzip_bits(resid):
    """Model-free compressed size (bits) of the integer residual sequence."""
    r = np.rint(resid).astype(np.int64)
    return 8.0 * len(gzip.compress(r.tobytes(), 9))


def emp_entropy_bits(resid):
    """K · H(empirical value distribution) -- marginal redundancy only."""
    r = np.rint(resid).astype(np.int64)
    if len(r) == 0:
        return 0.0
    _, counts = np.unique(r, return_counts=True)
    p = counts / counts.sum()
    H = -np.sum(p * np.log2(p))
    return float(H * len(r))


def _hinge_basis(u):
    """Linear-spline (hinge) basis with knots at every integer of u, plus extra
    knots in u∈[2,4] where Buchstab ω(u) bends most. O(max u)=O(log x) columns =
    polylog params. This captures ω(u)'s kink-at-every-integer structure that a
    low-degree polynomial cannot -- the FAIR polylog smooth model for phi(x,a)."""
    umin, umax = float(np.min(u)), float(np.max(u))
    knots = list(range(int(math.floor(umin)), int(math.ceil(umax)) + 1))
    knots += [2.5, 3.5]                       # extra resolution near the ω bend
    knots = sorted(set(k for k in knots if umin < k < umax))
    cols = [np.ones_like(u), u] + [np.maximum(0.0, u - k) for k in knots]
    return np.column_stack(cols), len(knots) + 2


def poly_fit_residual(phi, pls, x, degree):
    """Most-generous FAIR polylog predictor: fit phi/x as a linear spline in the
    Buchstab variable u=ln x/ln p_a with knots at integer u (O(log x) params --
    polylog). Captures the universal u·ω(u) Buchstab curve (kinked at integers);
    a global polynomial (the `degree` arg, kept for the ablation) cannot, leaving
    Θ(x/log x) smooth-model error. Returns (integer residual phi - x*fit, fit, npar)."""
    K = len(pls)
    u = math.log(x) / np.log(pls.astype(np.float64))      # u_a in [2, ln x/ln 2]
    y = phi[1:].astype(np.float64) / x                    # density in [0,1]
    B, npar = _hinge_basis(u)
    if K > npar + 1:
        coef, *_ = np.linalg.lstsq(B, y, rcond=None)
        fit = B @ coef
    else:                                                 # too few points: poly fallback
        deg = min(degree, K - 1)
        fit = np.polyval(np.polyfit(u, y, deg), u) if deg >= 1 else np.full(K, y.mean())
        npar = deg + 1
    pred = fit * x
    return phi[1:] - pred, pred, npar


def fit_exponent(xs, ys):
    """Least-squares slope of log2(y) vs log2(x); also vs √x reference."""
    lx = np.log2(np.asarray(xs, dtype=float))
    ly = np.log2(np.asarray(ys, dtype=float))
    A = np.vstack([lx, np.ones_like(lx)]).T
    slope, inter = np.linalg.lstsq(A, ly, rcond=None)[0]
    return float(slope), float(inter)


def build_residual_matrix(X0, spf, primes, window, degree):
    """Build the integer residual matrix R[i,a]=round(phi(X0+i,a) - smooth_fit) over
    a DENSE window i=0..window-1 of x, columns a=1..K (K=π(√(X0+window))). Incremental
    O(W·K); per-row polylog smooth fit subtracted (one shared spline basis at the
    window-midpoint u, since ln x is ~constant across the window -> multi-RHS lstsq).
    Returns (R, K, window, rrms, pls).

    Incremental: phi(x+1,a)=phi(x,a)+[a < idx(spf(x+1))] (x+1 is p_a-rough iff
    spf(x+1)>p_a)."""
    sqrtx = int(math.isqrt(X0 + window))
    K = int(np.searchsorted(primes, sqrtx, side='right'))
    pls = primes[:K]
    # base = phi(X0, a) for a=1..K (computed over THIS K, not phi_checkpoints' own)
    layer0 = np.searchsorted(primes, spf[2:X0 + 1])
    kill0 = np.bincount(layer0, minlength=len(primes))[:K]
    base = (X0 - np.cumsum(kill0)).astype(np.float64)  # length K
    # incremental build of M[i,a]=phi(X0+i,a): +1 to layers a<idx(spf(X0+i))
    idxs = np.searchsorted(primes, spf[X0 + 1: X0 + window])   # layer of spf for each new n
    aidx = np.arange(K)
    inc = (aidx[None, :] < idxs[:, None]).astype(np.float64)   # (window-1, K)
    M = np.empty((window, K), dtype=np.float64)
    M[0] = base
    M[1:] = base + np.cumsum(inc, axis=0)
    # subtract per-row polylog smooth fit; ln x ~constant across the window, so one
    # shared spline basis (midpoint u) fits all rows at once (multi-RHS lstsq).
    xcol = (X0 + np.arange(window)).astype(np.float64)[:, None]
    lnp = np.log(pls.astype(np.float64))
    uu = math.log(X0 + window // 2) / lnp
    B, _ = _hinge_basis(uu)
    Y = (M / xcol).T                                   # (K, window) densities
    coef, *_ = np.linalg.lstsq(B, Y, rcond=None)       # (npar, window)
    fit = (B @ coef).T                                 # (window, K) densities
    R = np.rint(M - fit * xcol)
    rrms = float(np.sqrt(np.mean(R ** 2)))
    return R, K, window, rrms, pls


def min_exact_rank(R):
    """Minimal SVD rank whose truncated reconstruction rounds to the EXACT integer
    matrix R (R already integer-valued, rint'd). Binary search on the spectrum.
    SVD subsumes the BEST low-rank (= best smooth) model, so this is the smooth-
    model-free integer degrees-of-freedom of R. rank ~ #cols => integer-incompressible
    (cert info-forced); rank polylog => compressible."""
    if R.size == 0:
        return 0
    U, S, Vt = np.linalg.svd(R, full_matrices=False)
    n = len(S)

    def exact_at(r):
        if r == 0:
            return bool(np.max(np.abs(R)) < 0.5)
        approx = (U[:, :r] * S[:r]) @ Vt[:r]
        return bool(np.max(np.abs(np.rint(approx) - R)) < 0.5)

    if not exact_at(n):
        return n                                   # not exact even at full numeric rank
    lo, hi = 0, n                                  # binary-search minimal exact rank
    while lo < hi:
        mid = (lo + hi) // 2
        if exact_at(mid):
            hi = mid
        else:
            lo = mid + 1
    return lo


def integer_rank(X0, spf, primes, window, degree):
    """Cross-check vs S65/CLOSED_PATHS row 737: minimal exact-integer-reconstruction
    rank of the smooth-removed checkpoint residual matrix over a dense window of x.
    Returns (rank, K, window, rrms). rank ~ K => integer-incompressible (cert
    info-forced); rank polylog => not."""
    R, K, window, rrms, pls = build_residual_matrix(X0, spf, primes, window, degree)
    return min_exact_rank(R), K, window, rrms


def layer_density_profile(R, K, nbins=10, fracs=(0.1, 0.25, 0.5, 0.75, 1.0)):
    """Is the joint √x information SPREAD across all K=π(√x) layers, or carried by
    o(K) of them? Two density statistics on the residual matrix R (W×K):

      * per-layer HARD-BITS profile, binned into nbins equal layer-index bins:
        mean ceil(log2(1+|r|)) per layer in each bin, plus the 'active fraction'
        (layers with >=1 hard bit) and the residual-ENERGY fraction. Uniform bits
        (every bin carries ~the same bits/layer) ⇒ dense; energy is scale-loaded
        (phi~x at small a) so bits/rank are the fair density measures.
      * PREFIX integer-reconstruction rank rank(R[:,:j]) at j=⌈f·K⌉: a LINEAR curve
        rank(j)≈(rank_K/K)·j means each block of layers adds fresh independent
        directions (dense); saturation (rank(K/2)≈rank(K)) means the rank — hence
        the info — lives in a vanishing prefix.

    Two scalar density metrics:
      bits_uniformity = min-bin / max-bin bits-per-layer   (→1 dense, →0 concentrated)
      rank_half_ratio = rank(R[:,:⌈K/2⌉]) / rank(R)         (→0.5 dense, →1 concentrated)
    """
    bits = np.where(np.abs(R) < 0.5, 0.0, np.ceil(np.log2(1.0 + np.abs(R))))
    colbits = bits.mean(axis=0)                       # mean bits/layer, per layer a
    colE = np.sum(R ** 2, axis=0)                     # residual energy, per layer a
    bins = np.array_split(np.arange(K), nbins)
    bin_bits = np.array([colbits[b].mean() for b in bins])
    bin_active = np.array([float(np.mean(colbits[b] > 0.5)) for b in bins])
    Etot = max(colE.sum(), 1e-30)
    bin_Efrac = np.array([colE[b].sum() / Etot for b in bins])
    rank_K = min_exact_rank(R)
    prefix = []
    for f in fracs:
        j = max(2, int(round(f * K)))
        prefix.append((f, j, min_exact_rank(R[:, :j])))
    rank_half = next((r for f, j, r in prefix if abs(f - 0.5) < 1e-9), rank_K)
    return dict(bin_bits=bin_bits, bin_active=bin_active, bin_Efrac=bin_Efrac,
                prefix=prefix, rank_K=rank_K,
                bits_uniformity=float(bin_bits.min() / max(bin_bits.max(), 1e-30)),
                rank_half_ratio=float(rank_half / max(rank_K, 1)),
                active_frac=float(np.mean(colbits > 0.5)))


# ----------------------------------------------------------------------------
# Main measurement over a dyadic sweep
# ----------------------------------------------------------------------------
def run_sweep(xmax, ks=None, degree=8, verbose=True, wfactors=(64, 192)):
    if ks is None:
        kmax = int(math.floor(math.log2(xmax)))
        ks = list(range(14, kmax + 1))
    if verbose:
        print(f"sieving smallest-prime-factor up to {xmax} ...", flush=True)
    spf, primes = sieve_spf(xmax)
    mu_small = mobius(40)

    rows = []
    for k in ks:
        x = 2 ** k
        if x > xmax:
            continue
        phi, K, sqrtx, pls = phi_checkpoints(x, spf, primes)
        pix = K + int(phi[K]) - 1                  # Legendre: π(x)=φ(x,K)+π(√x)-1
        # ---- predictors and their integer residuals over the K checkpoints ----
        mert = mertens_product(pls) * x            # (M) prime-aware
        asy = x * math.exp(-GAMMA) / np.log(pls.astype(np.float64))  # (A) polylog
        r_M = phi[1:] - mert
        r_A = phi[1:] - asy
        r_F, predF, nparF = poly_fit_residual(phi, pls, x, degree)   # (F) polylog fit
        # ---- reference upper bound (independent bits) ----
        H_M, H_A, H_F = magbits(r_M), magbits(r_A), magbits(r_F)
        # ---- DISCRIMINATOR: is the smooth-removed residual COMPRESSIBLE? ----
        gz_F = gzip_bits(r_F)                       # local/statistical redundancy
        ar_F = ar_residual_bits(r_F, order=3)       # neighbour-predictability
        ent_F = emp_entropy_bits(r_F)
        # per-layer bits profile: are surprises in few (late) layers or all of them?
        bits_layer = np.where(np.abs(np.rint(r_F)) < 0.5, 0.0,
                              np.ceil(np.log2(1.0 + np.abs(np.rint(r_F)))))
        half = len(bits_layer) // 2
        H_lo = float(bits_layer[:half].sum())       # first-half layers (small a)
        H_hi = float(bits_layer[half:].sum())       # second-half layers (large a)
        nz = int(np.count_nonzero(bits_layer))      # layers carrying >=1 hard bit
        # ---- single-value control (S36): π(x)-Li(x), π(x)-R(x) hard bits ----
        Lix = li(x) - li(2.0)
        Rx = riemann_R(x, mu_small)
        d_li = pix - Lix
        d_R = pix - Rx
        bits_li = math.log2(1 + abs(d_li))
        # ---- references ----
        sqrtx_bits = math.sqrt(x)                  # √x reference curve
        trivial = K * math.log2(2 + x)             # √x·log x trivial-cert proxy
        rows.append(dict(k=k, x=x, K=K, pix=pix, sqrtx=sqrtx,
                         H_M=H_M, H_A=H_A, H_F=H_F, gz_F=gz_F, ar_F=ar_F,
                         ent_F=ent_F, H_lo=H_lo, H_hi=H_hi, nz=nz,
                         nparF=nparF, d_li=d_li, d_R=d_R, bits_li=bits_li,
                         rF_rms=float(np.sqrt(np.mean(r_F**2))),
                         rF_max=float(np.max(np.abs(r_F))),
                         trivial=trivial, sqrtx_ref=sqrtx_bits))
        if verbose:
            print(f"  k={k:2d} x=2^{k}={x:>11d}  K=π(√x)={K:>5d}  π(x)={pix:>9d}"
                  f"  |π-Li|={abs(d_li):8.1f} ({bits_li:4.1f} bits)", flush=True)

    if verbose:
        _report(rows, degree)
        rank_sweep(spf, primes, xmax, degree, wfactors=wfactors)
    return rows


def _rank_sweep_pts(spf, primes, xmax, degree, wfactor, kmin=16, verbose=True):
    """Integer-reconstruction rank over a dyadic sweep of x, adaptive window
    W=wfactor*K. Returns the list of (x, K, rank). Window MUST scale with K
    (rank<=min(W,K); resolving all K modes needs W>>K)."""
    if verbose:
        print(f"\n=== integer-rank sweep (adaptive window W={wfactor}*K; S65/row-737 link) ===")
        print(f"{'k':>3} {'x':>12} {'K=π(√x)':>8} {'W':>8} {'rank':>6} {'rank/K':>7} {'rRMS':>9}")
    pts = []
    kmax = int(math.floor(math.log2(xmax)))
    for k in range(kmin, kmax + 1):
        X0 = 2 ** k
        K0 = int(np.searchsorted(primes, int(math.isqrt(X0)), side='right'))
        W = wfactor * K0
        if X0 + W >= xmax:
            continue
        rank, K, w, rrms = integer_rank(X0, spf, primes, W, degree)
        pts.append((X0, K, rank))
        if verbose:
            print(f"{k:>3} {X0:>12} {K:>8} {W:>8} {rank:>6} {rank/K:>7.2f} {rrms:>9.1f}")
    return pts


def _report_sweep_exponents(pts, wfactor):
    """α_rank, α_K, their ratio, and the rank/K floor over a sweep. The robust
    Θ(√x) statement is the rank/K FLOOR (rank>=c·K, no decay), independent of the
    finite-window log-log discount that depresses BOTH α_rank and α_K below 0.5."""
    if len(pts) < 3:
        print("  (too few points for an exponent fit)")
        return None
    xs = [p[0] for p in pts]
    s_rank, _ = fit_exponent(xs, [p[2] for p in pts])
    s_K, _ = fit_exponent(xs, [p[1] for p in pts])
    ratios = [p[2] / p[1] for p in pts]
    ratio_aa = s_rank / s_K if s_K else float('nan')
    print(f"  W={wfactor}*K:  α_rank={s_rank:+.3f}  α_K={s_K:+.3f}  "
          f"α_rank/α_K={ratio_aa:.3f}   rank/K∈[{min(ratios):.2f},{max(ratios):.2f}]")
    return dict(wfactor=wfactor, s_rank=s_rank, s_K=s_K, ratio_aa=ratio_aa,
                rk_min=min(ratios), rk_max=max(ratios), pts=pts)


def rank_sweep(spf, primes, xmax, degree, wfactors=(64, 192)):
    """Headline ROBUST measurement: integer-reconstruction rank of the smooth-removed
    checkpoint residual matrix over a dyadic sweep of x, at one or more adaptive window
    factors W=wfactor*K. SVD subsumes the BEST low-rank (= best smooth) model, so this
    is smooth-model-free.

    THE SHARPENING (this cycle): the finite-window α_rank tracks α_K (NOT 0.5) because
    K=π(√x)~2√x/ln x carries the SAME 1/ln x log-log discount -- so the honest Θ(√x)
    signal is (i) rank/K bounded below by a constant with NO downward drift, hence
    rank=Θ(K)=Θ(√x), and (ii) α_rank/α_K -> 1. A wider window (W=192K vs 64K) lifts the
    rank/K floor and removes the across-window rank/K climb that inflates α_rank at 64K,
    giving the cleaner α_rank≈α_K reading."""
    # window-sensitivity at one x: the convincing demonstration that the residual
    # is full-rank once adequately sampled (narrow windows under-count).
    kw = max(18, int(math.floor(math.log2(xmax))) - 4)
    X0w = 2 ** kw
    K0w = int(np.searchsorted(primes, int(math.isqrt(X0w)), side='right'))
    print(f"\n=== window-sensitivity of integer-rank at x=2^{kw} (K={K0w}) ===")
    print("  rank is W-capped until W >> K; full rank ⇒ K independent integer pieces:")
    for W in [K0w, 4 * K0w, 16 * K0w, 64 * K0w, 192 * K0w]:
        if X0w + W >= xmax:
            continue
        rank, K, w, rrms = integer_rank(X0w, spf, primes, W, degree)
        print(f"    W={W:>7} (W/K={W//K:>3})  rank={rank:>4}  K={K}  rank/K={rank/K:.2f}")

    summaries, pts_for_profile = [], None
    for wf in wfactors:
        pts = _rank_sweep_pts(spf, primes, xmax, degree, wf)
        s = _report_sweep_exponents(pts, wf)
        if s is not None:
            summaries.append(s)
        if pts:
            pts_for_profile = pts                      # keep the last (widest) sweep's pts

    # --- α_K reference over an EXTENDED dyadic range (no rank): shows BOTH α_rank and
    #     α_K -> 0.5 only as the PNT 1/ln x discount fades; the finite-window ~0.4 is
    #     the shared discount, not a sub-√x rank deficiency. ---
    kmax = int(math.floor(math.log2(xmax)))
    Kref = [(2 ** k, int(np.searchsorted(primes, int(math.isqrt(2 ** k)), side='right')))
            for k in range(16, kmax + 1)]
    if len(Kref) >= 3:
        sKfull, _ = fit_exponent([x for x, _ in Kref], [K for _, K in Kref])
        print(f"\n  α_K over full k=16..{kmax} = {sKfull:+.3f}  "
              f"(→0.5 as x→∞: d log K/d log x = 0.5 − 1/ln x + o(1); the finite-window")
        print("   reading is the PNT discount shared by rank — rank/K floor, not the")
        print("   exponent, is the sampling-robust √x statement)")

    print("\n  INTERPRETATION:")
    if summaries:
        best = summaries[-1]                           # widest window = most-sampled
        if best['rk_min'] > 0.80 and best['ratio_aa'] > 0.85:
            print(f"   Well-sampled rank ≈ {best['rk_min']:.2f}..{best['rk_max']:.2f} × K with NO")
            print(f"   downward drift, and α_rank/α_K={best['ratio_aa']:.2f}≈1: the K checkpoint")
            print("   residuals are integer-INDEPENDENT (no smaller smooth+low-rank model")
            print("   reconstructs them exactly), so rank=Θ(K)=Θ(√x)·polylog ⇒ the √x")
            print("   certificate is INFORMATION-FORCED for any sieve-reconstructing verifier")
            print("   -- a barrier, not mere construction shape.")
        elif best['s_rank'] > 0.30:
            print(f"   rank ~ x^{best['s_rank']:.2f}: super-polylog (no polylog cert) but possibly sub-√x.")
        else:
            print(f"   rank ~ x^{best['s_rank']:.2f}: (near-)polylog ⇒ √x cert NOT info-forced.")

    # --- per-layer DENSITY profile at the top of the (widest) sweep ---
    if pts_for_profile:
        X0p = pts_for_profile[-1][0]
        wf = wfactors[-1]
        K0p = int(np.searchsorted(primes, int(math.isqrt(X0p)), side='right'))
        Wp = wf * K0p
        if X0p + Wp < xmax:
            print(f"\n=== per-layer hard-bit DENSITY profile at x={X0p}=2^{int(round(math.log2(X0p)))} "
                  f"(W={Wp}, the widest sweep top) ===")
            R, K, w, rrms, pls = build_residual_matrix(X0p, spf, primes, Wp, degree)
            prof = layer_density_profile(R, K)
            print("  decile (by layer index a/K):  bits/layer  active-frac  energy-frac")
            for i, (bb, ba, be) in enumerate(zip(prof['bin_bits'], prof['bin_active'],
                                                 prof['bin_Efrac'])):
                bar = '#' * int(round(20 * bb / max(prof['bin_bits'].max(), 1e-9)))
                print(f"    dec {i}: {bb:6.2f}    {ba:5.2f}      {be:5.3f}  {bar}")
            print(f"  prefix integer-rank rank(R[:,:⌈fK⌉]):")
            for f, j, r in prof['prefix']:
                print(f"    f={f:.2f}  j={j:>4}  rank={r:>4}  rank/j={r / j:.2f}")
            print(f"\n  DENSITY metrics:")
            print(f"   bits_uniformity (min/max decile bits/layer) = {prof['bits_uniformity']:.2f}"
                  f"   (→1 uniform/dense; →0 concentrated)")
            print(f"   active fraction (layers with ≥1 hard bit)    = {prof['active_frac']:.2f}")
            print(f"   rank_half_ratio = rank(first K/2)/rank(K)    = {prof['rank_half_ratio']:.2f}"
                  f"   (→0.5 dense/linear; →1 carried by the prefix)")
            dense = (prof['bits_uniformity'] > 0.4 and prof['active_frac'] > 0.9
                     and 0.35 < prof['rank_half_ratio'] < 0.65)
            if dense:
                print("   ⇒ DENSE: the √x information is spread ~uniformly across all K=π(√x)")
                print("     layers (prefix rank grows LINEARLY, every layer carries ~equal hard")
                print("     bits), NOT carried by an o(K) subset -- so the joint info is genuinely")
                print("     Θ(K)=Θ(√x), not a handful of fat layers.")
            else:
                print("   ⇒ NON-uniform: inspect -- the √x may be carried by o(K) layers.")
    return pts_for_profile


def _report(rows, degree):
    xs = [r['x'] for r in rows]
    print("\n=== per-checkpoint residual (smooth-removed) and its compressibility ===")
    print(f"{'k':>3} {'K=π(√x)':>8} {'H_F=Σbits':>10} {'gzip_F':>9} {'AR3_F':>9} "
          f"{'rF_rms':>9} {'rF_max':>9} {'nz/K':>7} {'Hlo/Hhi':>9} {'√x':>9}")
    for r in rows:
        print(f"{r['k']:>3} {r['K']:>8} {r['H_F']:>10.0f} {r['gz_F']:>9.0f} "
              f"{r['ar_F']:>9.0f} {r['rF_rms']:>9.1f} {r['rF_max']:>9.0f} "
              f"{r['nz']/r['K']:>7.2f} {r['H_lo']:>4.0f}/{r['H_hi']:<4.0f} "
              f"{r['sqrtx_ref']:>9.0f}")

    print("\n=== fitted exponents  (slope of log2(quantity) vs log2(x);  √x ⇒ 0.5) ===")
    specs = [('H_F', 'Σ bits, predictor F (independent UPPER bound)', 'ref'),
             ('gz_F', 'gzip(residual_F)        -- DISCRIMINATOR', 'head'),
             ('ar_F', 'AR(3)-residual bits      -- DISCRIMINATOR', 'head'),
             ('ent_F', 'empirical-entropy bits (residual_F)', 'ref'),
             ('H_M', 'Σ bits, predictor M (prime-aware Mertens)', 'ref'),
             ('K', 'K=π(√x)                  reference (true 0.5)', 'ref'),
             ('rF_rms', 'residual_F RMS magnitude', 'ref'),
             ('trivial', 'trivial cert √x·log x    reference', 'ref')]
    exps = {}
    for key, label, _ in specs:
        ys = [r[key] for r in rows]
        if min(ys) <= 0:
            print(f"  {label:<46}: (nonpositive, skipped)")
            continue
        slope, _ = fit_exponent(xs, ys)
        exps[key] = slope
        print(f"  {label:<46}: alpha = {slope:+.3f}")

    print("\n=== single-value control (S36 / SESSION_INSIGHTS(i) reproduction) ===")
    bl = [r['bits_li'] for r in rows]
    slope_bl, _ = fit_exponent(xs, [max(b, 1e-9) for b in bl])
    print(f"  bits(|π(x)-Li(x)|): {bl[0]:.1f} .. {bl[-1]:.1f}  "
          f"(slope vs log2 x = {slope_bl:+.3f}; magnitude ~√x but BITS = O(log x) ≈ ½log2 x)")
    for r in rows:
        print(f"    x=2^{r['k']:<2}  π-Li={r['d_li']:+9.1f}  π-R={r['d_R']:+9.1f}"
              f"  ½log2x={0.5*r['k']:.1f}  bits={r['bits_li']:.1f}")

    slope_gz = exps.get('gz_F', 0.0)
    slope_ar = exps.get('ar_F', 0.0)
    slope_K = exps.get('K', 0.5)
    print("\n=== VERDICT ===")
    print(f"  K=π(√x) exponent          = {slope_K:+.3f}  (the cert has K blocks ⇒ SIZE is √x regardless;")
    print(f"                                       over a finite dyadic window √x reads as ~{slope_K:.2f}, not 0.50)")
    print(f"  Σ-bits (independent) exp  = {exps.get('H_F',0):+.3f}  (≈K-exp near-trivially: K layers × ≥1 bit)")
    print(f"  gzip(residual) exponent   = {slope_gz:+.3f}   <-- INFORMATION discriminator (info-forced ⟺ ≈ K-exp)")
    print(f"  AR(3)-residual exponent   = {slope_ar:+.3f}   <-- INFORMATION discriminator")
    disc = max(slope_gz, slope_ar)
    print("  (per-x is a single K-length sequence -- a 1D under-sample; the")
    print("   AUTHORITATIVE joint measure is the cross-x integer-rank sweep below.)")
    if disc > 0.6 * slope_K:
        print("\n  The smooth-removed residual is INCOMPRESSIBLE and scales like √x.")
        print("  ==> the chain-class √x certificate is INFORMATION-FORCED: the K layer")
        print("      checkpoints are ~K independent integer surprises (gzip & AR cannot")
        print("      shrink them), so NO certificate that reconstructs the sieve state can")
        print("      be sub-√x on information grounds. (Extends S65/row-737 from the φ")
        print("      matrix's integer-rank to the per-x checkpoint sequence's entropy.)")
    elif disc < 0.5 * slope_K:
        print("\n  The smooth-removed residual is COMPRESSIBLE (≪√x).")
        print("  ==> the √x certificate is NOT information-forced; the √x is the")
        print("      construction's sequential K-layer structure (S510), and a sub-√x")
        print("      certificate is not ruled out by incompressibility -- item-3's lower")
        print("      bound must be COMPUTATIONAL (#P-hardness), not information-theoretic.")
    else:
        print("\n  Intermediate compressibility -- inspect the per-layer profile & rank.")
    print("\n  (Distinct from the COMPUTATION barrier: a SINGLE π(x) value is O(log x)")
    print("   bits (control above, S36) yet √x-hard to COMPUTE; the CERTIFICATE question")
    print("   is the JOINT info of all K checkpoints, measured here.)")


# ----------------------------------------------------------------------------
# Selftest
# ----------------------------------------------------------------------------
def selftest():
    print("=== selftest ===")
    fails = 0

    def check(name, cond):
        nonlocal fails
        print(f"  [{'OK' if cond else 'FAIL'}] {name}")
        if not cond:
            fails += 1

    # 1. spf sieve correctness on small n
    spf, primes = sieve_spf(50)
    check("spf(2..10) correct",
          [int(spf[i]) for i in range(2, 11)] == [2, 3, 2, 5, 2, 7, 2, 3, 2])
    check("primes<=50 list", list(primes[:5]) == [2, 3, 5, 7, 11])

    # 2. phi boundary identities on a concrete x
    spf, primes = sieve_spf(2000)
    x = 1000
    phi, K, sqrtx, pls = phi_checkpoints(x, spf, primes)
    check("phi(x,0)=x", int(phi[0]) == x)
    check("phi(x,1)=x-floor(x/2)", int(phi[1]) == x - x // 2)
    # Legendre: π(x) = φ(x,K) + π(√x) - 1
    sieve_pi = int(np.searchsorted(primes, x, side='right'))
    pi_legendre = int(phi[K]) + K - 1
    check(f"Legendre π(1000)={pi_legendre}==sieve {sieve_pi}", pi_legendre == sieve_pi)

    # 3. phi matches a direct independent Legendre count (n=1 plus rough numbers)
    def phi_direct(x, a, pls):
        cnt = 0
        for n in range(1, x + 1):
            ok = True
            for j in range(a):
                if n % int(pls[j]) == 0:
                    ok = False
                    break
            if ok:
                cnt += 1
        return cnt
    ok_all = all(int(phi[a]) == phi_direct(x, a, pls) for a in [0, 1, 2, 3, K])
    check("phi(x,a) == direct Legendre count (a in {0,1,2,3,K})", ok_all)

    # 4. monotone non-increasing
    check("phi(x,a) non-increasing in a", bool(np.all(np.diff(phi) <= 0)))

    # 5. li/R sanity: Li(x) close to π(x), better than x/ln x; R even closer
    spf2, primes2 = sieve_spf(100000)
    x2 = 100000
    pix2 = int(np.searchsorted(primes2, x2, side='right'))   # π(10^5)=9592
    Lix2 = li(x2) - li(2.0)
    mu_small = mobius(40)
    Rx2 = riemann_R(x2, mu_small)
    check(f"π(10^5)=9592 (sieve gives {pix2})", pix2 == 9592)
    check(f"|π-Li|={abs(pix2-Lix2):.0f} < |π-x/lnx|",
          abs(pix2 - Lix2) < abs(pix2 - x2 / math.log(x2)))
    check(f"|π-R|={abs(pix2-Rx2):.0f} <= |π-Li|", abs(pix2 - Rx2) <= abs(pix2 - Lix2) + 1)

    # 6. metric estimators: constant seq -> ~0 bits; random seq -> ~K*entropy
    const = np.zeros(500)
    check("magbits(0-seq)=0", magbits(const) == 0.0)
    rng = np.random.default_rng(0)
    rnd = rng.integers(-100, 101, size=4000).astype(float)
    eb = emp_entropy_bits(rnd)
    # ~ uniform over ~201 values => ~log2(201)≈7.65 bits/symbol
    check(f"emp_entropy(random)≈K·7.6 (got {eb/4000:.2f} bits/sym)",
          7.0 < eb / 4000 < 8.2)
    gzc = gzip_bits(const)
    gzr = gzip_bits(rnd)
    check("gzip(const) << gzip(random)", gzc < 0.2 * gzr)

    # 7. exponent fitter recovers known slopes
    xs = [2.0 ** k for k in range(10, 25)]
    s_sqrt, _ = fit_exponent(xs, [math.sqrt(v) for v in xs])
    s_log, _ = fit_exponent(xs, [math.log2(v) for v in xs])
    check(f"fit_exponent(√x)≈0.5 (got {s_sqrt:.3f})", abs(s_sqrt - 0.5) < 1e-6)
    check(f"fit_exponent(log x)≈0 over range (got {s_log:.3f})", s_log < 0.10)

    # 8. poly_fit_residual: params polylog (<<K), residual integer-comparable
    r_F, pred, npar = poly_fit_residual(phi, pls, x, 8)
    check(f"fit uses polylog params npar={npar} << K={K}", npar < K)
    check("fit residual length == K", len(r_F) == K)

    # 9. The DISCRIMINATOR (gzip / AR) must SEPARATE incompressible from
    #    compressible K-length sequences -- magbits cannot (Σ over K layers ~√x
    #    for BOTH). Build, per x, an incompressible (iid) and a compressible
    #    (smooth deterministic) residual of length K, and check gzip exponents.
    xs2 = [2 ** k for k in range(14, 23)]
    Ks = [int(np.searchsorted(sieve_spf(int(math.isqrt(v)) + 1)[1],
                              int(math.isqrt(v)), side='right')) for v in xs2]
    rng = np.random.default_rng(1)
    gz_indep, gz_smooth, mb_indep, mb_smooth = [], [], [], []
    for v, Kv in zip(xs2, Ks):
        scale = math.sqrt(v)
        # INCOMPRESSIBLE: every one of K layers an independent wide surprise.
        indep = rng.integers(-int(scale), int(scale) + 1, size=Kv).astype(float)
        # COMPRESSIBLE: only polylog-many layers carry a surprise, the rest exact
        # (residual=0) -- the "sub-√x cert not ruled out" regime.
        smooth = np.zeros(Kv)
        nnz = max(2, int(math.log2(v)))
        idx = rng.choice(Kv, size=min(nnz, Kv), replace=False)
        smooth[idx] = rng.integers(-int(scale), int(scale) + 1, size=len(idx))
        gz_indep.append(gzip_bits(indep)); gz_smooth.append(gzip_bits(smooth))
        mb_indep.append(magbits(indep)); mb_smooth.append(magbits(smooth))
    sK, _ = fit_exponent(xs2, Ks)               # K=π(√x) exponent over this window
    si, _ = fit_exponent(xs2, gz_indep)
    ss, _ = fit_exponent(xs2, gz_smooth)
    mi, _ = fit_exponent(xs2, mb_indep)
    ms, _ = fit_exponent(xs2, mb_smooth)
    # info-forced (dense white): discriminators track K's exponent; compressible
    # (sparse): discriminators ~0. Thresholds are RELATIVE to sK (~0.42 here).
    sep = 0.6 * sK   # separator: info-forced (dense white) above, compressible below
    check(f"K=π(√x) exponent over window sK={sK:.2f}", 0.35 < sK < 0.50)
    check(f"gzip separates at {sep:.2f}: dense-white {si:.2f} > sep > sparse {ss:.2f}",
          si > sep > ss)
    check(f"magbits separates at {sep:.2f}: dense-white {mi:.2f} > sep > sparse {ms:.2f}",
          mi > sep > ms)

    # 10. integer_rank cross-check: synthetic exact-low-rank integer matrix has
    #     small rank; a full-entropy integer matrix needs near-full rank.
    A = np.outer(np.arange(40) % 7, np.arange(25) % 5).astype(float)  # rank<=2 exact
    U, S, Vt = np.linalg.svd(np.rint(A), full_matrices=False)
    r2 = next(r for r in range(1, len(S) + 1)
              if np.max(np.abs(np.rint((U[:, :r] * S[:r]) @ Vt[:r]) - A)) < 0.5)
    check(f"exact rank-2 integer matrix recovered at rank {r2}<=2", r2 <= 2)
    Rnd = rng.integers(-50, 51, size=(40, 25)).astype(float)
    U, S, Vt = np.linalg.svd(Rnd, full_matrices=False)
    rr = next((r for r in range(1, len(S) + 1)
               if np.max(np.abs(np.rint((U[:, :r] * S[:r]) @ Vt[:r]) - Rnd)) < 0.5),
              len(S))
    check(f"full-entropy integer matrix needs near-full rank {rr}>=20 (of 25)", rr >= 20)

    # 11. min_exact_rank REFACTOR SAFETY: the extracted helper reproduces the inline
    #     binary search, and integer_rank == build_residual_matrix + min_exact_rank.
    spf3, primes3 = sieve_spf(1 << 19)
    X0t = 1 << 16
    K0t = int(np.searchsorted(primes3, int(math.isqrt(X0t)), side='right'))
    Wt = 64 * K0t
    Rt, Kt, wt, rrmst, plst = build_residual_matrix(X0t, spf3, primes3, Wt, 8)
    rk_helper = min_exact_rank(Rt)
    rk_intfn, _, _, _ = integer_rank(X0t, spf3, primes3, Wt, 8)
    check(f"integer_rank == build+min_exact_rank (both {rk_helper})", rk_helper == rk_intfn)
    # min_exact_rank reproduces the inline-SVD binary search on an exact low-rank int matrix
    Alr = np.rint(np.outer(np.arange(60) % 9, np.arange(30) % 6)).astype(float)  # rank<=2
    U, S, Vt = np.linalg.svd(Alr, full_matrices=False)
    inline = next((r for r in range(0, len(S) + 1)
                   if np.max(np.abs(np.rint((U[:, :r] * S[:r]) @ Vt[:r]) - Alr)) < 0.5),
                  len(S))
    check(f"min_exact_rank matches inline on rank-2 matrix ({min_exact_rank(Alr)}=={inline})",
          min_exact_rank(Alr) == inline)
    check("min_exact_rank(zeros)=0", min_exact_rank(np.zeros((10, 5))) == 0)

    # 12. layer_density_profile SEPARATION: it must distinguish a DENSE residual
    #     (info spread across all K layers) from two kinds of CONCENTRATED residual.
    #     rank_half_ratio and (bits_uniformity, active_frac) are COMPLEMENTARY
    #     discriminators -- the conjunction is what the verdict uses.
    rng2 = np.random.default_rng(7)
    Wd, Kd = 2000, 80
    # (a) DENSE: all K columns independent integer surprises.
    R_dense = rng2.integers(-300, 301, size=(Wd, Kd)).astype(float)
    pd = layer_density_profile(R_dense, Kd)
    check(f"dense: rank_half≈0.5 (got {pd['rank_half_ratio']:.2f})",
          0.35 < pd['rank_half_ratio'] < 0.65)
    check(f"dense: active≈1 ({pd['active_frac']:.2f}) & bits uniform ({pd['bits_uniformity']:.2f}>0.4)",
          pd['active_frac'] > 0.9 and pd['bits_uniformity'] > 0.4)
    # (b) RANK-concentrated: second half = exact COPIES of first-half columns. Every
    #     column is "active" with uniform bits, but the rank lives entirely in the
    #     first K/2 -> rank_half_ratio≈1. Only the rank metric catches this.
    half = Kd // 2
    R_rc = np.empty((Wd, Kd))
    R_rc[:, :half] = rng2.integers(-300, 301, size=(Wd, half)).astype(float)
    R_rc[:, half:] = R_rc[:, :half]                      # duplicate -> no new rank
    prc = layer_density_profile(R_rc, Kd)
    check(f"rank-concentrated caught: rank_half>0.9 (got {prc['rank_half_ratio']:.2f}) "
          f"though bits uniform ({prc['bits_uniformity']:.2f})",
          prc['rank_half_ratio'] > 0.9)
    # (c) BIT-concentrated: only ~log2(K) columns carry surprises, the rest exact (0).
    R_bc = np.zeros((Wd, Kd))
    nnz = max(2, int(math.log2(Kd)))
    cols = rng2.choice(Kd, size=nnz, replace=False)
    R_bc[:, cols] = rng2.integers(-300, 301, size=(Wd, nnz)).astype(float)
    pbc = layer_density_profile(R_bc, Kd)
    check(f"bit-concentrated caught: active<0.5 ({pbc['active_frac']:.2f}) & "
          f"bits_uniformity<0.1 ({pbc['bits_uniformity']:.2f})",
          pbc['active_frac'] < 0.5 and pbc['bits_uniformity'] < 0.1)

    print(f"\n{'ALL PASS' if fails == 0 else str(fails)+' FAILED'}")
    return fails == 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--xmax', type=int, default=1 << 22,
                    help='sieve bound / largest x=2^k in the sweep (default 2^22)')
    ap.add_argument('--degree', type=int, default=8,
                    help='polylog smooth-fit polynomial degree (predictor F)')
    ap.add_argument('--kmin', type=int, default=14)
    ap.add_argument('--wfactors', type=int, nargs='+', default=[64, 192],
                    help='adaptive rank-window factors W=wfactor*K to sweep '
                         '(default 64 192: the 64->192 sharpening of α_rank/α_K)')
    args = ap.parse_args()

    if args.selftest:
        ok = selftest()
        sys.exit(0 if ok else 1)

    kmax = int(math.floor(math.log2(args.xmax)))
    ks = list(range(args.kmin, kmax + 1))
    run_sweep(args.xmax, ks=ks, degree=args.degree, wfactors=tuple(args.wfactors))


if __name__ == '__main__':
    main()
