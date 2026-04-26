"""
Fourth-encoding search (FOCUS-6).

Enumerate additive / multiplicative number-theoretic functions f: N -> Q whose
summatory S_f(x) = sum_{n<=x} f(n) might be informationally equivalent to pi(x)
up to polylog conversion.

For each candidate, we measure:
  (a) Polylog computability of S_f(x): does it admit a known closed form
      (Dirichlet hyperbola / known Dirichlet series factoring through zeta) so
      that S_f(x) is evaluable in O(polylog(x)) time? -- ANALYTIC, hand-coded.
  (b) Empirical residual R_f(x) = S_f(x) - M_f(x) where M_f(x) is the analytic
      main term, computed by least-squares fit to the canonical asymptotic
      basis when the closed form is not pre-coded.
  (c) Growth rate alpha s.t. |R_f| ~ x^alpha   (from log-log slope on a sample
      grid in [1e3, 1e5]).
  (d) Pearson correlation rho(R_f, E_pi)  where  E_pi(x) = pi(x) - Li(x).
      |rho| -> 1  ==>  R_f reduces to the explicit-formula error (mode E).
      |rho| -> 0  ==>  no prime information beyond the smooth part (mode I).
  (e) Free-identity probes: is S_f(x) mod p trivial (constant in x or
      = poly(x) mod p) for p in {2, 3, 5}? This catches Liouville-style
      L(x) mod 2 = x mod 2 traps before they look promising.
  (f) Empirical classification into one of {C, E, I, NOVEL}.

The candidates were selected to span: additive, multiplicative, totient-like,
divisor-like, smooth-/rough-counting, partition-flavoured, digit-flavoured,
and analytic-special-function-flavoured number-theoretic functions.

Range: x in {1, 2, ..., X_MAX = 100_000}.  Computation is dominated by an
O(X_MAX log X_MAX)-style sieve of small prime factors and a single linear scan
to accumulate every f(n).  Total wall time on a laptop: < 30 s.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass

import numpy as np
import sympy
from mpmath import li, mp


X_MAX = 100_000
PROBE_MIN = 1_000  # for log-log slope and rho computation
PROBE_GRID = np.unique(
    np.round(np.geomspace(PROBE_MIN, X_MAX, 240)).astype(np.int64)
)


# ---------------------------------------------------------------------------
# 1.  Sieve of arithmetic primitives we reuse across candidates.
# ---------------------------------------------------------------------------

def build_primitives(N: int):
    """Build per-n vectors of standard arithmetic functions for n=0..N.

    Returns a dict with: lpf, mu, lam, omega_big, Omega_big, sigma0, sigma1,
                          phi, is_prime, smallest_prime_pow_exp, J2.
    """
    print(f"  building primitives sieve up to N={N} ...", flush=True)
    t0 = time.time()

    lpf = np.zeros(N + 1, dtype=np.int64)        # least prime factor
    Omega_big = np.zeros(N + 1, dtype=np.int64)  # with multiplicity
    omega_big = np.zeros(N + 1, dtype=np.int64)  # distinct
    mu = np.ones(N + 1, dtype=np.int64)
    sigma0 = np.zeros(N + 1, dtype=np.int64)
    sigma1 = np.zeros(N + 1, dtype=np.int64)
    phi = np.arange(N + 1, dtype=np.int64)
    J2 = np.arange(N + 1, dtype=np.int64) ** 2

    # build lpf, mu via standard linear sieve
    is_prime = np.ones(N + 1, dtype=bool)
    is_prime[:2] = False
    for p in range(2, N + 1):
        if is_prime[p]:
            for m in range(p, N + 1, p):
                if lpf[m] == 0:
                    lpf[m] = p
                if m != p:
                    is_prime[m] = False

    # mu, omega_big, Omega_big via factorisation using lpf
    for n in range(2, N + 1):
        m = n
        squarefree = True
        omega = 0
        Omega = 0
        last_p = 0
        while m > 1:
            p = lpf[m]
            count = 0
            while m % p == 0:
                m //= p
                count += 1
            omega += 1
            Omega += count
            if count >= 2:
                squarefree = False
            last_p = p
        mu[n] = (-1) ** omega if squarefree else 0
        omega_big[n] = omega
        Omega_big[n] = Omega
    mu[1] = 1
    omega_big[1] = 0
    Omega_big[1] = 0

    # sigma0, sigma1 directly from prime factorisation (cheap enough at N=1e5)
    for n in range(1, N + 1):
        m = n
        s0 = 1
        s1 = 1
        ph = 1
        j2 = 1
        while m > 1:
            p = lpf[m]
            e = 0
            while m % p == 0:
                m //= p
                e += 1
            s0 *= (e + 1)
            s1 *= (p ** (e + 1) - 1) // (p - 1)
            ph *= p ** (e - 1) * (p - 1)
            j2 *= p ** (2 * (e - 1)) * (p * p - 1)
        sigma0[n] = s0
        sigma1[n] = s1
        phi[n] = ph
        J2[n] = j2

    sigma0[0] = 0
    sigma1[0] = 0
    phi[0] = 0
    J2[0] = 0

    lam = ((-1) ** Omega_big).astype(np.int64)  # Liouville
    lam[0] = 0

    print(f"  primitives built in {time.time()-t0:.1f}s.", flush=True)
    return {
        "lpf": lpf,
        "mu": mu,
        "lam": lam,
        "omega_big": omega_big,
        "Omega_big": Omega_big,
        "sigma0": sigma0,
        "sigma1": sigma1,
        "phi": phi,
        "J2": J2,
        "is_prime": is_prime,
    }


# ---------------------------------------------------------------------------
# 2.  Candidate definitions.
#
# Each candidate provides:
#   name                : short identifier
#   per_n               : np.ndarray (N+1,) giving f(n)
#   polylog_evaluable   : bool, "is sum_{n<=x} f(n) known to be O(polylog)?"
#   smooth_form         : str, the analytic main term (for the report)
#   notes               : str
# ---------------------------------------------------------------------------

def make_candidates(prims, N: int) -> list[dict]:
    n = np.arange(N + 1, dtype=np.float64)

    chi_P = prims["is_prime"].astype(np.int64)

    # log Gamma fractional part: f(n) = {log Gamma(n+1)} (fractional part).
    # Sum is x log x - x + ... + sum of fractional residues.  Per-n: log(n)
    # is the natural smooth part; we use the fractional residual.
    log_n = np.log(np.where(n == 0, 1.0, n))
    # cumulative integer log: floor(log Gamma(n+1)) is unwieldy; instead use
    # the fluctuating part log(n) - mean(log).  This isolates the
    # fluctuating contribution.
    log_n_resid = log_n.copy()
    log_n_resid[0] = 0.0

    # 1/n harmonic
    inv_n = np.zeros(N + 1, dtype=np.float64)
    inv_n[1:] = 1.0 / n[1:]

    # smooth-number indicator at B=20: 1 if all prime factors <= 20
    B_smooth = 20
    smooth_indicator_B20 = np.zeros(N + 1, dtype=np.int64)
    smooth_indicator_B20[0] = 0
    smooth_indicator_B20[1] = 1
    for k in range(2, N + 1):
        m = k
        ok = True
        while m > 1:
            p = prims["lpf"][m]
            if p > B_smooth:
                ok = False
                break
            while m % p == 0:
                m //= p
        smooth_indicator_B20[k] = 1 if ok else 0

    # rough-number indicator at B=20: 1 if smallest prime factor > 20
    rough_indicator_B20 = ((prims["lpf"] > B_smooth) | (n == 1)).astype(np.int64)
    rough_indicator_B20[0] = 0

    # digit sum base 10
    s10 = np.zeros(N + 1, dtype=np.int64)
    for k in range(1, N + 1):
        x = k
        s = 0
        while x:
            s += x % 10
            x //= 10
        s10[k] = s

    # binary digit sum (Hamming weight)
    s2 = np.zeros(N + 1, dtype=np.int64)
    for k in range(1, N + 1):
        s2[k] = bin(k).count("1")

    # p-adic valuation v_2(n)
    v2 = np.zeros(N + 1, dtype=np.int64)
    for k in range(1, N + 1):
        v2[k] = (k & -k).bit_length() - 1

    # r_2(n) = #ways to write n = a^2 + b^2 with a,b in Z (signed, ordered)
    r2 = np.zeros(N + 1, dtype=np.int64)
    for a in range(-int(math.isqrt(N)) - 1, int(math.isqrt(N)) + 2):
        a2 = a * a
        if a2 > N:
            continue
        for b in range(-int(math.isqrt(N - a2)) - 1, int(math.isqrt(N - a2)) + 2):
            s = a2 + b * b
            if 0 <= s <= N:
                r2[s] += 1

    # largest prime factor LPF(n)
    LPF = np.zeros(N + 1, dtype=np.int64)
    for k in range(2, N + 1):
        m = k
        last = 0
        while m > 1:
            p = prims["lpf"][m]
            last = p
            while m % p == 0:
                m //= p
        LPF[k] = last

    # smallest prime factor minus 1 (additive on prime-power components only crudely)
    lpf_minus_one = np.where(prims["lpf"] > 0, prims["lpf"] - 1, 0)

    # exp Stieltjes-like fluctuation: f(n) = log(p_n) - log(n) for n prime, 0 else
    # this is *prime-dependent*, kept as a control (mode C diagnostic).

    # Lambda (von Mangoldt): log p if n=p^k, 0 else
    Lambda = np.zeros(N + 1, dtype=np.float64)
    for p in range(2, N + 1):
        if prims["is_prime"][p]:
            pk = p
            while pk <= N:
                Lambda[pk] = math.log(p)
                pk *= p

    # Mertens-like signed: lam(n)/n  (Liouville weighted)
    lam_over_n = np.zeros(N + 1, dtype=np.float64)
    lam_over_n[1:] = prims["lam"][1:] / n[1:]

    candidates = [
        {
            "name": "chi_P (prime indicator)",
            "per_n": chi_P.astype(np.float64),
            "polylog": False,
            "smooth_form": "Li(x) (~ x/log x)",
            "notes": "trivial control: S_chi_P(x) = pi(x). Cost = mode C/E/I depending on method.",
        },
        {
            "name": "Lambda (von Mangoldt)",
            "per_n": Lambda,
            "polylog": False,
            "smooth_form": "x",
            "notes": "psi(x). Reduces to explicit formula; same complexity as pi(x).",
        },
        {
            "name": "lambda (Liouville)",
            "per_n": prims["lam"].astype(np.float64),
            "polylog": False,
            "smooth_form": "0  (sublinear, |L(x)| << x^{1/2+eps})",
            "notes": "Closed S55 (free identity L mod 2 = x mod 2). Includes for control.",
        },
        {
            "name": "mu (Mobius)",
            "per_n": prims["mu"].astype(np.float64),
            "polylog": False,
            "smooth_form": "0",
            "notes": "Mertens M(x). Helfgott-Thompson O(x^{3/5}); transfer to pi closed S16.",
        },
        {
            "name": "Omega (factor count w/ mult.)",
            "per_n": prims["Omega_big"].astype(np.float64),
            "polylog": False,
            "smooth_form": "x log log x + B*x + ...",
            "notes": "Erdos-Kac. Sum has Mertens-type error.",
        },
        {
            "name": "omega (distinct factor count)",
            "per_n": prims["omega_big"].astype(np.float64),
            "polylog": False,
            "smooth_form": "x log log x + B*x",
            "notes": "Sum tied to PNT remainder term (Selberg).",
        },
        {
            "name": "sigma_0 (#divisors)",
            "per_n": prims["sigma0"].astype(np.float64),
            "polylog": False,
            "smooth_form": "x log x + (2*gamma-1)*x + Delta(x)",
            "notes": "Dirichlet divisor problem. Delta(x) is its OWN open problem, ~x^{1/4}.",
        },
        {
            "name": "sigma_1 (sum of divisors)",
            "per_n": prims["sigma1"].astype(np.float64),
            "polylog": False,
            "smooth_form": "(pi^2/12) * x^2 + O(x log x)",
            "notes": "Sum dominated by smooth quadratic term.",
        },
        {
            "name": "phi (Euler totient)",
            "per_n": prims["phi"].astype(np.float64),
            "polylog": False,
            "smooth_form": "(3/pi^2) * x^2 + O(x log x)",
            "notes": "Mertens 1874. Error term tied to Mertens function.",
        },
        {
            "name": "J_2 (Jordan totient)",
            "per_n": prims["J2"].astype(np.float64),
            "polylog": False,
            "smooth_form": "x^3 / (3 * zeta(3))",
            "notes": "Higher-order analogue; Dirichlet series = zeta(s-2)/zeta(s).",
        },
        {
            "name": "log n",
            "per_n": log_n_resid,
            "polylog": True,
            "smooth_form": "x log x - x + (1/2)log(2*pi*x) + 1/(12x) - ... (Stirling)",
            "notes": "Stirling asymptotic; closed form in O(polylog).",
        },
        {
            "name": "1/n (harmonic)",
            "per_n": inv_n,
            "polylog": True,
            "smooth_form": "log x + gamma + 1/(2x) - 1/(12 x^2) + ...",
            "notes": "Hurwitz zeta; closed form in O(polylog).",
        },
        {
            "name": "20-smooth indicator",
            "per_n": smooth_indicator_B20.astype(np.float64),
            "polylog": True,
            "smooth_form": "rho-function * x  (Dickman; B fixed)",
            "notes": "B=20 fixed: Psi(x,20) computable in polylog. But B=poly(log x) for any prime info.",
        },
        {
            "name": "20-rough indicator",
            "per_n": rough_indicator_B20.astype(np.float64),
            "polylog": False,
            "smooth_form": "x * prod_{p<=20}(1-1/p)",
            "notes": "Wheel sieve. B=sqrt(x) recovers pi via Legendre identity.",
        },
        {
            "name": "digit sum base 10",
            "per_n": s10.astype(np.float64),
            "polylog": True,
            "smooth_form": "(9/2) * x * log_10(x) + ...",
            "notes": "Trollope-Delange; closed form. Independent of primes (mode I).",
        },
        {
            "name": "binary digit sum (popcount)",
            "per_n": s2.astype(np.float64),
            "polylog": True,
            "smooth_form": "(1/2) * x * log_2(x) + x*F(log_2 x) (Trollope)",
            "notes": "Trollope. F is bounded fractal; independent of primes (mode I).",
        },
        {
            "name": "v_2 (2-adic valuation)",
            "per_n": v2.astype(np.float64),
            "polylog": True,
            "smooth_form": "x - s_2(x)",
            "notes": "sum_{n<=x} v_2(n) = x - s_2(x). Independent of primes (mode I).",
        },
        {
            "name": "r_2 (#repr a^2+b^2)",
            "per_n": r2.astype(np.float64),
            "polylog": False,
            "smooth_form": "pi * x  (Gauss circle problem)",
            "notes": "Sum = #lattice points in disk of radius sqrt(x). Error O(x^{1/3}).",
        },
        {
            "name": "LPF (largest prime factor)",
            "per_n": LPF.astype(np.float64),
            "polylog": False,
            "smooth_form": "x^2 / (2 log x) (Erdos-Tenenbaum)",
            "notes": "Sum involves Dickman; computing pointwise = factoring.",
        },
        {
            "name": "lpf - 1",
            "per_n": lpf_minus_one.astype(np.float64),
            "polylog": False,
            "smooth_form": "x^2 / (2 (log x)^2) (rough)",
            "notes": "Smallest-prime-factor based; mode C (factor lookup).",
        },
        {
            "name": "lambda(n)/n",
            "per_n": lam_over_n,
            "polylog": False,
            "smooth_form": "0",
            "notes": "Liouville weighted; D-series = zeta(2s)/zeta(s).",
        },
    ]

    return candidates


# ---------------------------------------------------------------------------
# 3.  Empirical analysis pipeline.
# ---------------------------------------------------------------------------

def fit_smooth_basis(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Least-squares fit y(x) ~ c0 + c1*x + c2*x*log(x) + c3*log(x) + c4*x^2.

    Returns (residual, coeffs).
    """
    X = np.column_stack(
        [
            np.ones_like(x, dtype=np.float64),
            x.astype(np.float64),
            x.astype(np.float64) * np.log(x.astype(np.float64) + 1e-12),
            np.log(x.astype(np.float64) + 1e-12),
            (x.astype(np.float64) ** 2),
        ]
    )
    coeffs, *_ = np.linalg.lstsq(X, y.astype(np.float64), rcond=None)
    pred = X @ coeffs
    return y - pred, coeffs


def loglog_growth_alpha(probe_x: np.ndarray, probe_R: np.ndarray) -> float:
    """Fit log|R| ~ alpha*log(x) + b on the largest 60% of probes."""
    mask = (probe_R != 0) & (probe_x > 1)
    if mask.sum() < 10:
        return float("nan")
    lo = int(mask.sum() * 0.4)
    sx = np.log(probe_x[mask])[lo:]
    sy = np.log(np.abs(probe_R[mask]))[lo:]
    A = np.column_stack([sx, np.ones_like(sx)])
    sol, *_ = np.linalg.lstsq(A, sy, rcond=None)
    return float(sol[0])


def free_identity_check(per_n: np.ndarray, S_f_int: np.ndarray, X: np.ndarray) -> dict:
    """Probe whether S_f(x) mod p is a polynomial in x mod p for p in {2,3,5}.

    Returns dict { 'mod_p_explained_R2': ... } where explained_R2 is the variance
    explained by a degree-2 polynomial in x, taken modulo p.
    """
    out: dict[str, float] = {}
    if not np.all(np.isfinite(per_n)):
        return {"mod_2": float("nan"), "mod_3": float("nan"), "mod_5": float("nan")}

    if not np.allclose(per_n, np.round(per_n)):
        return {"mod_2": float("nan"), "mod_3": float("nan"), "mod_5": float("nan")}

    for p in (2, 3, 5):
        Smod = S_f_int.astype(np.int64) % p  # over the probe grid in X
        S_at_X = Smod[X]
        x = X.astype(np.float64)
        # try to fit S mod p as round((a + b x + c x^2) mod p)
        # easier: just compute mutual reduction R = corr(S_at_X mod p, x mod p)
        x_modp = X % p
        # exact-equality fraction:
        out[f"mod_{p}"] = float((S_at_X == x_modp).mean())
    return out


@dataclass
class CandidateResult:
    name: str
    polylog: bool
    growth_alpha: float
    rho_with_E_pi: float
    R2_against_smooth: float
    free_id_mod2: float
    free_id_mod3: float
    free_id_mod5: float
    notes: str
    classification: str


def classify(growth_alpha: float, rho: float, polylog: bool, name: str) -> str:
    """Heuristic triage based on residual scaling and explicit-formula coupling."""
    if "chi_P" in name or "Lambda" in name:
        return "C/E (control)"
    abs_rho = abs(rho)
    if not polylog:
        # not polylog -> can't be a fourth informationally-complete primitive
        # but we still record whether the residual encodes prime info
        if abs_rho > 0.5:
            return "E (zeta-zero residual)"
        if growth_alpha < 0.4:
            return "I (smooth, info loss)"
        return "E? (prime-coupled, sub-sqrt)"
    # polylog: candidate-of-interest only if strongly correlated AND informative
    if abs_rho > 0.6 and growth_alpha > 0.4:
        return "*** CANDIDATE FOURTH ENCODING ***"
    if abs_rho < 0.1:
        return "I (polylog but pi-independent)"
    return "I (polylog, weak prime coupling)"


def main() -> int:
    mp.dps = 30
    print("== Fourth-encoding search ==")
    print(f"X_MAX = {X_MAX}, probe grid size = {len(PROBE_GRID)}")

    prims = build_primitives(X_MAX)
    candidates = make_candidates(prims, X_MAX)

    # canonical reference: pi(x) and Li(x) on the probe grid
    print("  computing pi(x) on full range ...", flush=True)
    pi_full = np.cumsum(prims["is_prime"].astype(np.int64))
    li_probe = np.array([float(li(int(x))) for x in PROBE_GRID])
    E_pi_probe = (pi_full[PROBE_GRID].astype(np.float64) - li_probe).astype(np.float64)

    results: list[CandidateResult] = []

    print()
    for cand in candidates:
        name = cand["name"]
        per_n = cand["per_n"]
        S_f_full = np.cumsum(per_n)

        # Fit smooth basis
        resid_full, coeffs = fit_smooth_basis(np.arange(X_MAX + 1), S_f_full)
        S_at_probe = S_f_full[PROBE_GRID]
        resid_probe = resid_full[PROBE_GRID]

        # variance explained by smooth basis
        var_resid = float(np.var(resid_full[PROBE_MIN:]))
        var_total = float(np.var(S_f_full[PROBE_MIN:]))
        R2_smooth = 1.0 - var_resid / max(var_total, 1e-30)

        # log-log growth slope of |residual|
        alpha = loglog_growth_alpha(PROBE_GRID, resid_probe)

        # correlation with E_pi
        if np.std(resid_probe) > 0 and np.std(E_pi_probe) > 0:
            rho = float(np.corrcoef(resid_probe, E_pi_probe)[0, 1])
        else:
            rho = 0.0

        # free identity check (only for integer-valued candidates)
        per_n_int = np.allclose(per_n, np.round(per_n))
        if per_n_int:
            S_int = np.cumsum(per_n.astype(np.int64))
            free = free_identity_check(per_n, S_int, PROBE_GRID)
        else:
            free = {"mod_2": float("nan"), "mod_3": float("nan"), "mod_5": float("nan")}

        cls = classify(alpha, rho, cand["polylog"], name)

        results.append(
            CandidateResult(
                name=name,
                polylog=cand["polylog"],
                growth_alpha=alpha,
                rho_with_E_pi=rho,
                R2_against_smooth=R2_smooth,
                free_id_mod2=free["mod_2"],
                free_id_mod3=free["mod_3"],
                free_id_mod5=free["mod_5"],
                notes=cand["notes"],
                classification=cls,
            )
        )

        print(
            f"  {name:38s}  polylog={cand['polylog']!s:5s}  "
            f"alpha={alpha:+.3f}  rho={rho:+.4f}  "
            f"R2_smooth={R2_smooth:.6f}  "
            f"free(mod2,3,5)=({free['mod_2']:.3f},{free['mod_3']:.3f},{free['mod_5']:.3f})  "
            f"-> {cls}"
        )

    # ---- summary ----
    print()
    novel = [r for r in results if "CANDIDATE" in r.classification]
    print(f"Novel-fourth-encoding hits: {len(novel)}")
    for r in novel:
        print(f"  *** {r.name}: rho={r.rho_with_E_pi:+.4f}, alpha={r.growth_alpha:+.3f}")

    # save machine-readable summary
    out_path = (
        "/apps/aplikacijos/prime-research/experiments/algebraic/"
        "fourth_encoding_search/fourth_encoding_search_data.csv"
    )
    with open(out_path, "w") as f:
        f.write(
            "name,polylog,growth_alpha,rho_with_E_pi,R2_against_smooth,"
            "free_id_mod2,free_id_mod3,free_id_mod5,classification\n"
        )
        for r in results:
            f.write(
                f"\"{r.name}\",{r.polylog},{r.growth_alpha:+.4f},"
                f"{r.rho_with_E_pi:+.6f},{r.R2_against_smooth:.6f},"
                f"{r.free_id_mod2:.4f},{r.free_id_mod3:.4f},{r.free_id_mod5:.4f},"
                f"\"{r.classification}\"\n"
            )
    print(f"\nWrote {out_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
