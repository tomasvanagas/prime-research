"""
Session 4: ZERO SUM COMPRESSION

The KEY observation from all explorations:
  π(x) = li(x) - Σ_ρ li(x^ρ) + small terms

The barrier is summing over ALL zeros ρ = 1/2 + iγ_k.

BUT: The sum S(x) = Σ_k li(x^{1/2+iγ_k}) is a DETERMINISTIC function of x.
It's not random - it just LOOKS random because the γ_k are incommensurate.

QUESTION: Can we find a CLOSED FORM or FAST ALGORITHM for S(x)?

Key ideas:
A) S(x) = Σ_k f(γ_k, x) where f is known. If Σ_k f(γ_k, x) can be evaluated
   as a contour integral involving ζ'/ζ, we avoid enumerating zeros.

B) Actually, the explicit formula IS this contour integral:
   π(x) = (1/2πi) ∫_{c-i∞}^{c+i∞} log ζ(s) * x^s/s ds
   This integral CAN be evaluated... but it still costs O(x^{1/2+ε}).

C) What about APPROXIMATING S(x) via the argument principle?
   N(T) = #{γ ≤ T} = T/(2π) ln(T/(2πe)) + 7/8 + S(T)
   where S(T) = (1/π) arg ζ(1/2 + iT) (the S-function)

D) BACKLUND'S EXACT FORMULA: The Riemann-von Mangoldt formula gives
   N(T) EXACTLY as a smooth term + arg ζ. Can we use a similar trick
   for the sum Σ li(x^ρ)?

E) PARTIAL SUMMATION on the zero sum:
   Σ_{γ≤T} li(x^{1/2+iγ}) = N(T)*li(x^{1/2+iT}) - ∫_0^T N(t) d/dt[li(x^{1/2+it})] dt

   N(t) is SMOOTH (modulo S(t) fluctuations). So the integral is computable!
   But... d/dt li(x^{1/2+it}) = x^{1/2+it} * i*ln(x) / ((1/2+it)*ln(x))
   This oscillates wildly, making the integral hard to evaluate.

F) STATIONARY PHASE: The integral ∫ N(t) * x^{it} * g(t) dt
   has stationary phase when d/dt[t*ln(x)] = 0, i.e., never (no stationary point).
   So the integral oscillates and cancels... which means the TAIL of the zero
   sum (γ > T) should be small. But we already know that - the error is x/T.

Let me try a completely different angle:

G) DIGIT EXTRACTION via the BBP formula analogy.
   The Bailey-Borwein-Plouffe formula computes the k-th hex digit of π
   without computing preceding digits. This works because π has a special
   series representation.

   Question: Does π(x) or p(n) have a BBP-like representation?

   For π(x), the "digits" in some base could be related to...
   Modular arithmetic on the explicit formula.

H) THE NUCLEAR OPTION: What if we accept ~50% digits from R^{-1}(n)
   and try to get the REMAINING digits from a different, independent source?

   R^{-1}(n) gives the TOP half of digits of p(n).
   The BOTTOM half comes from the oscillatory correction.

   Can we compute the bottom half INDEPENDENTLY using:
   - p(n) mod small primes (residue information)
   - Density of primes in specific arithmetic progressions
   - Some other source of information about p(n)?
"""

import mpmath
import math
import time

mpmath.mp.dps = 50

# Let's test idea H: Can residue information help?
#
# We know from PNT in arithmetic progressions:
# π(x; q, a) ~ li(x) / φ(q) for (a,q)=1
#
# The ERROR in this approximation is ~x^{1/2} (under GRH).
# So knowing π(x; q, a) mod φ(q) is as hard as knowing π(x).
#
# BUT: what about HEURISTIC/STATISTICAL information?
# We know p(n) ≠ 0 mod p for any prime p (except p itself).
# This gives us p(n) mod 2 = 1 (for n > 1), p(n) mod 3 ∈ {1,2} (for n > 2), etc.
# Not useful - just says the candidate is coprime to small primes.

# Let's try idea A/D more carefully.
# The EXACT explicit formula (under RH):
# ψ(x) = x - Σ_ρ x^ρ/ρ - ln(2π) - (1/2)ln(1-x^{-2})  [for x not a prime power]
#
# Σ_ρ x^ρ/ρ = 2 Re Σ_{γ>0} x^{1/2+iγ}/(1/2+iγ)
#            = 2√x Re Σ_{γ>0} x^{iγ}/(1/2+iγ)
#            = 2√x Re Σ_{γ>0} e^{iγ ln x}/(1/2+iγ)
#
# Define F(u) = Σ_{γ>0} e^{iγu}/(1/2+iγ)  where u = ln x
#
# F(u) is a trigonometric series with "frequencies" γ_k (the zeta zeros).
#
# Can F(u) be computed WITHOUT enumerating γ_k?
#
# By the argument principle:
# F(u) = (1/2πi) ∮ e^{isu}/(1/2+is) * (ζ'/ζ)(1/2+is) ds ... not quite
#
# Actually: Σ_ρ f(ρ) = (1/2πi) ∮ f(s) * (-ζ'/ζ)(s) ds
# where the contour encloses all zeros.
#
# So: Σ_{γ>0} x^{1/2+iγ}/(1/2+iγ) = (1/2πi) ∮ x^s/s * (-ζ'/ζ)(s) ds
#                                     (contour enclosing zeros in upper half)
#
# This is EXACTLY the contour integral for ψ(x)!
# The contour integral for ψ(x) = (1/2πi) ∫_{c-i∞}^{c+i∞} (-ζ'/ζ)(s) x^s/s ds
#
# Evaluating this contour integral by residues gives the explicit formula.
# But evaluating it NUMERICALLY gives ψ(x) directly!
#
# The contour integral can be evaluated in O(T * polylog(T)) time
# where T is the integration range. For exact ψ(x), need T ~ x.
# But truncating at T ~ x^{1/2} gives error ~ x/T = x^{1/2}.
# This is the analytic method of Lagarias-Odlyzko!

# SO: the contour integral is NOT a shortcut. It's equivalent to the analytic method.

# Let me try yet another angle:
# IDEA I: Self-referential / Fixed-point approach
#
# Define the map T: x → li^{-1}(n + Σ_{ρ} li(x^ρ) + ln2 - ...)
# The fixed point of T is x = p(n) (the nth prime).
#
# Start with x_0 = R^{-1}(n).
# x_1 = T(x_0) involves evaluating Σ_ρ li(x_0^ρ) which requires ALL zeros.
# Unless... we truncate to K zeros and iterate.
#
# Each iteration with K zeros: the "wrong" zeros contribute error ~x/(γ_K)
# to the correction, but the ITERATION might amplify or dampen this.
#
# Key question: is T contracting? If |T'(x*)| < 1, iteration converges.
# T'(x) ≈ d/dx [li^{-1}(n + Σ li(x^ρ) + ...)]
#         = (1/li'(x*)) * Σ_ρ d/dx li(x^ρ)
#         = x*ln(x*) * Σ_ρ x^{ρ-1}/ln(x^ρ)
#         ≈ Σ_ρ x^{ρ-1} * x*ln(x)/(ρ*ln x)
#         = Σ_ρ x^ρ/(ρ*x) * (ln x)/(... )
# This is O(√x * K / x) = O(K/√x)
# For K < √x, |T'| < 1 and the iteration CONVERGES!
# But the fixed point is the WRONG one - it's the fixed point of the
# TRUNCATED map, not the full map.

# IDEA J: QUANTUM COMPUTATION CONNECTION
# Shor's algorithm factors N in O(polylog N) using quantum period-finding.
# If we could factor efficiently, we could build π(x) from factoring...
# But even quantum factoring gives factorizations, not prime counts.
#
# A quantum algorithm for π(x): embed π(x) as the eigenvalue of some operator.
# The "quantum counting" problem. Actually, Grover search gives O(√N) speedup.
# For π(x) via analytic method: quantum speedup gives O(x^{1/4+ε}). Still too slow.

# OK let me actually try to compute something novel.
# IDEA K: FAST PARTIAL ZERO SUM via the Z-function

# The Riemann-Siegel Z-function: Z(t) = e^{iθ(t)} ζ(1/2 + it)
# where θ(t) = arg Γ(1/4 + it/2) - (t/2) ln π
#
# Z(t) is real-valued and Z(t) = 0 iff ζ(1/2+it) = 0.
#
# The Riemann-Siegel formula gives Z(t) ≈ 2 Σ_{n≤√(t/2π)} n^{-1/2} cos(θ(t) - t ln n)
# This takes O(√t) operations.
#
# For the zero sum Σ_{γ≤T} f(γ), we can use:
# Σ_{γ≤T} f(γ) = ∫_0^T f(t) dN(t) where N(t) = (smooth) + S(t)/π
# = ∫_0^T f(t) (smooth part)' dt + ∫_0^T f(t) dS(t)/π
#
# The smooth part: ∫ f(t) * [ln(t/2π)/(2π)] dt - computable analytically!
# The oscillatory part: ∫ f(t) dS(t)/π - this is where the difficulty is.
#
# S(t) = (1/π) arg ζ(1/2+it) fluctuates like O(ln t).
# Can we bound or compute ∫ f(t) dS(t) efficiently?

# Actually, let me try a PRACTICAL experiment:
# Compute the zero sum using only the SMOOTH part of N(T) and see how good it is.

print("=" * 70)
print("ZERO SUM COMPRESSION: Can we avoid enumerating zeros?")
print("=" * 70)

# First, let's see what the "smooth zero sum" gives us.
# Replace Σ_{γ_k} f(γ_k) with ∫ f(t) * N'_smooth(t) dt
# where N_smooth(t) = t/(2π) * ln(t/(2πe)) + 7/8

def N_smooth(t):
    """Smooth part of zero counting function N(T)."""
    if t <= 0:
        return mpmath.mpf(0)
    t = mpmath.mpf(t)
    return t/(2*mpmath.pi) * mpmath.log(t/(2*mpmath.pi*mpmath.e)) + mpmath.mpf(7)/8

def N_smooth_deriv(t):
    """Derivative of smooth N(T)."""
    t = mpmath.mpf(t)
    return mpmath.log(t/(2*mpmath.pi)) / (2*mpmath.pi)

def zero_sum_smooth(x, T_max=1000):
    """
    Approximate Σ_{γ>0} li(x^{1/2+iγ})/(1/2+iγ) using smooth N(T).

    Σ f(γ) ≈ ∫_0^T f(t) N'(t) dt (smooth approximation)
    where f(t) = x^{1/2+it} / ((1/2+it) * ln(x) * (1/2+it))
    Actually f(t) = Ei((1/2+it) ln x)
    """
    x = mpmath.mpf(x)
    lnx = mpmath.log(x)

    def integrand(t):
        t = mpmath.mpf(t)
        if t < 0.1:
            return mpmath.mpf(0)
        rho = mpmath.mpc(0.5, t)
        # li(x^ρ) = Ei(ρ * ln x)
        val = mpmath.ei(rho * lnx)
        n_prime = N_smooth_deriv(t)
        return mpmath.re(val) * n_prime

    # Numerical integration
    result = mpmath.quad(integrand, [1, T_max], error=True)
    return 2 * result[0]  # factor 2 for conjugate pairs

def pi_smooth_zero_sum(x, T_max=1000):
    """π(x) using smooth zero sum approximation."""
    x = mpmath.mpf(x)
    result = mpmath.li(x) - mpmath.log(2)
    correction = zero_sum_smooth(x, T_max)
    result -= correction

    # Möbius inversion (just first few terms)
    # This is Π(x), need to subtract Π(√x)/2 etc.
    sqrtx = mpmath.sqrt(x)
    if sqrtx > 2:
        result_sq = mpmath.li(sqrtx) - mpmath.log(2) - zero_sum_smooth(sqrtx, T_max)
        result -= result_sq / 2

    cbrtx = mpmath.power(x, mpmath.mpf(1)/3)
    if cbrtx > 2:
        result_cb = mpmath.li(cbrtx) - mpmath.log(2) - zero_sum_smooth(cbrtx, T_max)
        result -= result_cb / 3

    return result

# Reference primes
def sieve(n):
    s = [True] * (n+1)
    s[0] = s[1] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]:
            for j in range(i*i, n+1, i):
                s[j] = False
    return [i for i, v in enumerate(s) if v]

primes = sieve(120000)

print("\n--- Smooth zero sum approximation of π(x) ---")
print("Replacing discrete sum over zeros with integral against smooth N'(T)")
print()

for n in [50, 100, 200, 500, 1000]:
    x = primes[n-1]
    t0 = time.time()
    pi_val = float(pi_smooth_zero_sum(x, T_max=200))
    elapsed = time.time() - t0
    error = abs(pi_val - n)
    print(f"  n={n:5d}, x={x:7d}: π_smooth(x)={pi_val:.4f}, "
          f"error={error:.4f}, time={elapsed:.2f}s")

# Compare with discrete zeros
ZEROS = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
         37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
         52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
         67.079811, 69.546402, 72.067158, 75.704691, 77.144840]

def pi_discrete_zeros(x, K=20):
    """π(x) using discrete sum over K known zeros."""
    x = mpmath.mpf(x)
    result = mpmath.li(x) - mpmath.log(2)
    lnx = mpmath.log(x)

    for gamma in ZEROS[:K]:
        gamma = mpmath.mpf(gamma)
        rho = mpmath.mpc(0.5, gamma)
        li_val = mpmath.ei(rho * lnx)
        result -= 2 * mpmath.re(li_val)

    # Möbius inversion
    sqrtx = mpmath.sqrt(x)
    if float(sqrtx) > 2:
        r2 = mpmath.li(sqrtx) - mpmath.log(2)
        for gamma in ZEROS[:K]:
            gamma = mpmath.mpf(gamma)
            rho = mpmath.mpc(0.5, gamma)
            r2 -= 2 * mpmath.re(mpmath.ei(rho * mpmath.log(sqrtx)))
        result -= r2 / 2

    return result

print("\n--- Comparison: smooth integral vs discrete zeros ---")
for n in [100, 500, 1000]:
    x = primes[n-1]
    smooth_val = float(pi_smooth_zero_sum(x, T_max=200))
    discrete_val = float(pi_discrete_zeros(x, K=20))
    print(f"  n={n:5d}: smooth={smooth_val:.4f} (err={abs(smooth_val-n):.4f}), "
          f"discrete={discrete_val:.4f} (err={abs(discrete_val-n):.4f})")

print("\n--- CONCLUSION ---")
print("The smooth zero sum replaces Σ f(γ_k) with ∫ f(t) N'(t) dt.")
print("If this works, it avoids enumerating zeros entirely!")
print("But N'(t) misses the oscillatory part S(t), causing errors.")
print("The oscillatory part S(t) = (1/π) arg ζ(1/2+it) is exactly")
print("what encodes the irregularity of primes - it cannot be avoided.")
