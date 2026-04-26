# Fresh-Perspective Brainstorm — Session 53 (2026-04-26 evening)

Spirit: *novel angles not in earlier brainstorms, not in the wildcard catalog,
not in `proposals/`*. Each item names the concrete falsifier.

A separate brainstorm earlier today
(`fresh_perspective_session_2026_04_26.md`) listed 10 angles; items 3–10
remain untested. This session picks one of those (Item 3, free-probability
moments) AND adds the seven below.

---

## A. Discrete Schrödinger ground-state count = π(x)?

Build the tridiagonal matrix `H[n,n] = log n`, `H[n,n±1] = -1` for
`n ∈ [1, x]`. By the Sturm sequence theorem the number of eigenvalues
below `E` is computable in `O(x)` from sign changes of the principal
minors. **Question:** does the count of eigenvalues below `E* = log x − 1`
(a natural "Fermi level") equal `π(x)` plus polylog correction?

If yes, Sturm sequencing gives a path to fast counting. If no, the
mismatch tells us how far purely-positional potentials are from
encoding primality.

**Falsifier:** compute eigenvalue counts at small x and compare to π(x).

## B. Tropical eigenvalue of `M[i,j] = log(p_{i+j})` matrix

In `(max, +)` algebra, the eigenvalue of a square matrix is the maximum
mean of any cycle (Karp). Build the `k × k` matrix where `M[i,j] = log
p_{i+j}` (primes indexed). Compute its tropical eigenvalue λ(k). If
λ(k) has a closed form in `k`, we get prime-position information out.

Most likely: λ(k) = log(p_{2k-1}) by the diagonal trace dominance, which
is just an obfuscated prime lookup. Worth checking.

**Falsifier:** compute λ(k) for k ≤ 20; check for non-trivial closed form.

## C. Sheaf Euler characteristic on Spec ℤ slice

Define the constant sheaf ℤ on the open subset `U_x = Spec ℤ ∖
{(p) : p > x}` of `Spec ℤ`. Its Euler characteristic
`χ(U_x, ℤ)` should encode the cardinality of the closed subscheme
{(p) : p ≤ x}, i.e., π(x). Question: does a derived-category
spectral sequence give χ in `O(polylog x)` operations?

This is a reformulation rather than a new algorithm. Unless the
spectral sequence collapses (which would be a miracle), the
computation reduces to enumerating primes.

**Falsifier:** read derived-category literature; find any non-trivial
collapse for `Spec ℤ` slice sheaves. Currently a thought experiment;
no clear empirical test.

## D. Lee–Yang zeros of an "Ising-with-Coulomb" prime model

Spin variable `σ_n = +1` if `n` is prime, `-1` else, on `n ∈ [1, x]`.
Energy `H = − Σ_{m<n} σ_m σ_n / (n − m)`. Partition function
`Z_x(β) = ⟨exp(−βH)⟩` over **uniform** spin configurations. The
fugacity-zeros (Lee–Yang) of `Z_x(β)` lie on the unit circle for
ferromagnetic Ising. If they lie on a sparse arc, **zero count = π(x)
counts a structural quantity that polylog reads from arc geometry**.

**Falsifier:** compute Z_x(β) on small x ≤ 16 by brute summation over
2^x configurations; locate Lee–Yang zeros numerically; check if their
density at the critical point relates to π(x).

## E. Multiplicative complexity of `1_P` over `GF(2^k)` for k > 1

It is known the AND-circuit complexity of "n is prime" over `GF(2)`
is super-polynomial (project's circuit-complexity findings). Over a
larger field `GF(2^k)`, AND gates compute a richer function in one
step. **Test:** count AND gates needed to implement `1_P` on n-bit
inputs as a circuit over `GF(2^k)` for `k = 2, 3, 4`. If the gate
count drops by a factor super-polynomial in `k`, we have evidence
for fast field-arithmetic primality testing.

**Falsifier:** small-N enumeration via MQ-encoding; existing MQ
solvers can find shortest circuits up to ≈ 10 inputs.

## F. Random projection mod random small primes

Pick `q_1, …, q_T` uniformly random primes ≤ Q. For each `q_i`,
compute `π(x) mod q_i` via *some other algorithm whose mod-q_i answer
is fast*. Combine via CRT.

The bottleneck is the sub-question: is there ANY family of primes
`q` such that `π(x) mod q` is computable in `O(polylog x)`?
The project's pseudorandomness findings on `π(x) mod 2` suggest
no, but the test for `q = 3, 5, 7, 11, 13` has not been done.

**Falsifier:** explicit DFA/automaticity tests on `π(x) mod q`.

## G. Frequency-comb representation of δ(x)

The explicit formula gives `δ(x) ~ −Σ_ρ x^ρ/ρ`. View this as a
"frequency comb" with frequencies `Im(ρ_k) / log x` and amplitudes
`x^{1/2}/|ρ_k|`. **Hypothesis:** the comb has only `O(polylog x)`
*non-degenerate* frequencies once we account for symmetries
(`ρ_k ↔ 1−ρ̄_k`) AND clustering from GUE level repulsion.

**Falsifier:** for `x = 10^k`, compute the truncated zero-sum at
truncation `T = 2^j` for `j = 5, 10, 15, 20`. Plot RMS error vs.
`T`. If RMS drops faster than `1/√T` (the GUE prediction), we have
non-trivial compression. If RMS ~ `1/√T`, no compression. (This
overlaps Item G in the earlier brainstorm but at a different scale.)

---

## This session's executed experiment

**Item 3 from earlier brainstorm: free-probability moments of
δ(x)/√x.** The free Wigner / free Poisson / Marchenko-Pastur
distributions have known moment sequences (Catalan-related). If the
empirical moments of normalized δ on dyadic windows match a free
distribution, the implied structure could yield an R-transform-based
shortcut.

Code: `free_probability_delta.py`.
