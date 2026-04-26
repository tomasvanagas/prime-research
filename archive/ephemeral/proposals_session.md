# Proposal Session — 2026-04-26 (Fresh, no prior-context)

Goal: O(polylog(n)) exact computation of p(n).
Constraint: do not look at CLAUDE.md / CLOSED_PATHS.md before generating.
Discipline: every proposal lists math idea, pseudocode, complexity, key
assumption, and an n < 10000 test.

---

## Proposal A — Tropical Hankel Rank of π

### Idea

Build the Hankel-like matrix
$$H \in \mathbb{N}^{m \times m}, \qquad H_{i,j} \;=\; \pi(i+j),\quad 0 \le i, j < m,$$
where $m = \lfloor N/2 \rfloor$. Treat $H$ as a matrix over the
**tropical (max-plus)** semiring $\mathbb R \cup \{-\infty\}$ with
$\oplus = \max$, $\otimes = +$.

The **tropical rank** $r(H)$ is the smallest $r$ such that $H = U \otimes V$
with $U \in (\mathbb{R}\cup\{-\infty\})^{m \times r}$,
$V \in (\mathbb{R}\cup\{-\infty\})^{r\times m}$. In contrast to ordinary linear
algebra, tropically low rank means $\pi(i+j) = \max_k (u_{ik} + v_{kj})$.
Such a decomposition is exactly a piecewise-linear concave-in-each-variable
representation — and the prime-counting function is asymptotically piecewise
linear (with breakpoints at primes).

If $r(H_N)$ scales like $\mathrm{polylog}(N)$, then primes can be located by
a sequence of $r$ comparisons of linear forms; combined with binary search,
this gives $O(\mathrm{polylog}\, N)$ evaluation of $\pi(x)$, hence $p(n)$
by inversion.

### Pseudocode

```python
def tropical_rank_pi(N):
    pi = sieve_of_eratosthenes_pi(2*N+2)        # baseline
    m  = N
    H  = np.array([[pi[i+j] for j in range(m)]
                              for i in range(m)], dtype=float)
    return cgb_tropical_rank(H)                 # see below
```

The **Cuninghame–Green–Butkovic** algorithm reduces $H$ by repeatedly
finding a column that is a max-plus combination of others; complexity is
$O(m^4)$ but only needed for the *test*, not the actual evaluator.

### Complexity

If $r(H_N) = O(\log^c N)$ then evaluating $\pi(x)$ for $x \le N$ takes
$O(r \log m)$ comparisons after $O(r m)$ preprocessing.
For $p(n)$: invert via tropical bisection in $O(r \log^2 N)$.

### Key assumption

Tropical rank of $\pi$-Hankel is sub-polynomial, not full.
Equivalent (heuristic): the prime gap landscape, after $\log$-lift, is
piecewise-linear with $\mathrm{polylog}$ many pieces.

### Test for n < 10000

* Build $H$ of size $m=2000$ ($N=4000$) and $m=5000$ ($N=10000$).
* Compute tropical rank via CGB.
* If the rank curve $r(N)$ fits $A \log^c N$ with $c \le 3$, advance.
* If $r(N)$ scales polynomially (e.g., $r \approx N^{1/3}$), CLOSED.

---

## Proposal B — Continued-Fraction Spectrum of the Prime Constant

### Idea

Define the **prime constant**
$$\alpha \;=\; \sum_{p\,\text{prime}} 2^{-p} \;=\; 0.414682509851111660248109622\ldots_{(2)}$$
i.e. $\alpha = 0.b_1 b_2 b_3 \ldots$ in binary with $b_k = 1$ iff $k$ is
prime. The classical real-line representation of $\alpha$ is
information-equivalent to $\pi(\cdot)$.

Switch representations: compute the **regular continued fraction**
$\alpha = [a_0; a_1, a_2, \ldots]$. The classical Lévy theorem says for
random reals $\frac{1}{n}\log a_n \to \pi^2 / (12 \log 2)$ (Khintchine
constant), but the prime constant is *not* random — it is computable.

**Conjecture (testable, fresh):** the partial-quotient sequence
$\{a_n\}$ admits a sub-linear sketch — e.g. its DFT decays like $1/k$ or
the sequence is automatic in a base related to small primes.

If true, we can compute $a_n$ in $\mathrm{polylog}\, n$ time, then
recover the binary digits (i.e., the primality indicator up to position
$\sim n$) by reconstructing convergents — yielding $\pi(x)$ in
$\mathrm{polylog}\, x$.

### Pseudocode

```python
def cf_of_prime_constant(N_digits):
    sieve   = sieve_of_eratosthenes(N_digits)
    mp.prec = 4 * N_digits
    alpha   = mp.mpf(0)
    for p in sieve:
        alpha += mp.mpf(2) ** (-p)
    a, x = [], alpha
    for _ in range(M):                          # M ≪ N_digits depth
        ai = int(mp.floor(x))
        a.append(ai)
        if x == ai:        break
        x = 1 / (x - ai)
    return a
```

Then probe `a` for structure:
* power-spectrum decay (FFT of `log(a)`),
* automaticity (Allouche–Shallit test),
* PSLQ on the generating series $A(z) = \sum a_n z^n$,
* k-automatic structure for small k.

### Complexity

If `a_n` is $k$-automatic, $a_n$ is computable in $O(\log n)$ (one DFA
walk) and convergents in $O(\log^2 n)$, giving polylog primality test.
Otherwise the experiment shows a barrier (Khintchine-like growth ⇒ random).

### Key assumption

Some non-Khintchine, non-random structure is hiding in the CF expansion
of $\alpha$. Most $\alpha$'s are Khintchine-typical; the prime
constant is special — it's deterministic.

### Test for n < 10000

* Compute CF of $\alpha$ from $\sim 12000$ binary digits → expect ~7000
  partial quotients (Khintchine-on-average gives ~1.7 digits per step).
* Apply: power-spectrum, autocorrelation, automaticity tests.
* CLOSE if $\{a_n\}$ shows full Khintchine statistics with no structural
  deviation.

---

## Proposal C — Algebraicity of Riemann Zero Differences

### Idea

The Riemann explicit formula needs $O(\sqrt{x})$ zeros to compute
$\pi(x)$ exactly. But what if the zeros themselves are not independent?

**Conjecture C (folklore-violating):** there exists a
$d = \mathrm{polylog}\,N$ -dimensional $\mathbb{Q}$-vector space
$V_d \subset \mathbb{R}$ such that $\{\gamma_k\}_{k \le N} \subset V_d$
modulo errors $< 2^{-N}$ (i.e., the zero ordinates lie on a
finite-rank lattice up to polynomially small noise).

If true, only $d$ "atomic" zero ordinates suffice to reconstruct all
$\sqrt{x}$ contributions, dropping the explicit-formula cost from
$\sqrt{x}$ to $d$.

This is the Hilbert–Pólya conjecture from the data-compression angle:
HP says zeros are eigenvalues of a Hermitian operator $H$. If $H$ has a
$\mathrm{polylog}$-dimensional algebraic shadow, $\pi$ is polylog.

### Pseudocode

```python
def zero_basis_search(zeros, basis):
    # zeros: γ_1,...,γ_N
    # basis: (1, π, log2, log3, log5, ..., log p_k,
    #        log Γ-values, multiple zeta values, ...)
    relations = []
    for k in range(1, len(zeros)):
        delta_k = zeros[k] - zeros[0]
        rel     = pslq([delta_k] + list(basis), tol=1e-40)
        if rel:                                  # found integer relation
            relations.append((k, rel))
    return relations
```

Then check:
* what fraction of $\gamma_k$ admit relations of bounded height,
* if all relations factor through a small sub-basis (then we have $V_d$).

### Complexity

Pre-cache $d$ atomic zeros, compute $\pi(x)$ via explicit formula in
$O(d \cdot \log^c x)$.

### Key assumption

The zeros are NOT $\mathbb{Q}$-linearly independent (folklore expects
they ARE). PSLQ would find any relation of polynomial height; absence
gives a strong negative result, presence is a breakthrough.

### Test for n < 10000

* Use the existing 8000 zeros at 30-digit precision.
* For $k = 2, \ldots, 200$, run PSLQ on $\gamma_k - \gamma_1$ against
  basis $\{1, \pi, \log 2, \log 3, \log 5, \log 7, \log 11, \log
  \Gamma(1/4), \log \Gamma(1/3)\}$.
* Record any candidate relation; bound height.
* If 0/200 relations found at height $\le 10^{15}$, proposal C is dead.
* If $\ge 5/200$ found, escalate immediately to proposal C-1 (full-basis
  search) and write to `novel/`.

---

## Proposal D — Reservoir Echo of π(x) via Liquid State Machines

### Idea

A **reservoir computer** is a fixed random recurrent net; only the
output linear layer is trained. Universal-approximation-on-trajectories
results (Maass–Natschläger–Markram 2002, Grigoryeva–Ortega 2018) say
that an $L$-node reservoir with sufficient memory and fading memory can
approximate any $\epsilon$-fading-memory functional of its input.

If $\pi$ has fading memory in $\log n$ (i.e., $\pi(n)$ is determined by
recent inputs only), an $L = \mathrm{polylog}\,n$ reservoir suffices.

If we *succeed* in fitting $L = O(\log^c n)$ reservoir to $\pi$ on $n
\le 10000$ with zero error and stable extrapolation to $10001\ldots
20000$, we have empirical evidence that $\pi$ has polylog descriptive
complexity.

### Pseudocode

```python
def reservoir_pi(N_train, N_test, L):
    rng = np.random.default_rng(0)
    W   = rng.standard_normal((L, L)) / np.sqrt(L)
    Win = rng.standard_normal((L, 1))
    rho = 0.9 / max(np.abs(np.linalg.eigvals(W)))
    W  *= rho
    state = np.zeros(L)
    Xs    = []
    for n in range(1, N_train + N_test + 1):
        u     = np.array([float(n)])
        state = np.tanh(W @ state + Win @ u)
        Xs.append(state.copy())
    Xtr, Xte = Xs[:N_train], Xs[N_train:]
    ytr      = pi_table[1:N_train+1]
    Wout, *_ = np.linalg.lstsq(np.array(Xtr), ytr, rcond=None)
    pred     = np.array(Xte) @ Wout
    return pred, pi_table[N_train+1:N_train+N_test+1]
```

### Complexity

Inference: $O(L^2 \log n)$ per query. If $L = \log^c n$, polylog.

### Key assumption

$\pi$ has fading memory (recent values determine current up to bounded
correction). Empirically falsifiable in one experiment.

### Test for n < 10000

* Train on $n \in [1, 5000]$, test on $n \in [5001, 10000]$.
* Sweep $L \in \{4, 8, 16, 32, 64, 128, 256, 512\}$.
* If extrapolation MAE stays $\ge 1$ across $L$, fading-memory
  assumption fails — close.
* If for some $L = O(\log^c 10000)$ MAE = 0, escalate.

---

## Triage / order of execution

1. **Proposal C** — cheapest (data exists), highest signal: yes/no on
   whether zeros admit hidden algebraic structure.
2. **Proposal B** — moderate cost, fresh angle, gives information either
   way (Khintchine vs. structure).
3. **Proposal A** — moderate cost, novel framing (tropical), buildable.
4. **Proposal D** — clearest empirical pass/fail, but most likely to
   degenerate into "deep nets memorise but don't generalise."

Run order: C → B → A → D.
