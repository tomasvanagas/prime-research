# Brandt 2024 MKtP technique vs pi(x) mod 2 — formal definitions

## Source

Nicholas Brandt, "Lower Bounds for Levin–Kolmogorov Complexity," TCC 2024,
LNCS 15363 (Springer), and IACR Cryptology ePrint Archive 2024/687.

## The objects

### Brandt's hypothesis class (the function class for which the lower bound applies)

The Brandt 2024 lower bound is **specifically** the language

```
    MKtP_U  :=  { (y, k) in {0,1}^m x [m]  :  Kt_U(y) <= k }
```

where `U` is a fixed prefix-free universal Turing machine (UTM) and

```
    Kt_U(x)  :=  min { |Pi| + ceil(log2 t)  :  U(Pi) outputs x in <= t steps }.
```

Brandt's main theorem (Theorem 1) is

> `MKtP_U  not in  DTIME[O(n)]`     (unconditional)

with corollaries:
- `MKtP_U  not in  DTIME[O(n^2)]` if the underlying machine model has a
  linear-time universal simulation (e.g. RAMs).
- `MKtP_U  not in  Heur_{0, gamma_fn} DTIME[O(n)]` for one-sided heuristics
  with small false-negative rate (Theorem 2).
- `MKtP_U  not in  Union_k DTIME[ prod_{i=0..k} ln^{(i)}(n) ]` (slightly
  super-linear, Corollary 2).

### Candidate target

The *natural* target this project cares about is `pi mod 2 : N -> {0,1}`.
For each N we form

```
    pi_N  :=  ( pi(0) mod 2,  pi(1) mod 2,  ...,  pi(2^N - 1) mod 2 ) in {0,1}^{2^N}.
```

This is a single string per N, computed deterministically from x by sieving.
The associated boolean function f_N : {0,1}^N -> {0,1} reads the input x
in binary and outputs pi(x) mod 2.

### The conjectured relationship (this experiment tests)

The hope: there exists a reduction or proof-skeleton transplant such that

> "any TM Pi_pi solving pi(x) mod 2 in linear (or sub-quadratic) time
>  yields a TM Pi_TRA(Pi_pi) that, by a Brandt-style traversal,
>  produces a string z with Kt(z) >= |z| using only the time budget of
>  Pi_pi, contradicting the Kt definition."

The four-clause structural argument in `brandt_mktp.py` shows this is
**false**: the traversal is welded to Kt itself (Brandt page 5: "we need
to exploit some property exhibited by the actual set R_Kt of Kt-random
strings but not by any set R_PiBB"), and 1-Kt-randomness of Chaitin's
Omega (Lemma 2) is the unique ingredient that bypasses the black-box
barrier. None of these ingredients appear in any natural-function
hypothesis.

## Signature of the bounded-Kt simulator used in `brandt_mktp.py`

We implement a 3-bit per-instruction stack VM with these ops:

| Bits | Mnemonic | Effect |
|------|----------|--------|
| 000  | PUSH 0   | emit '0' |
| 001  | PUSH 1   | emit '1' |
| 010  | DUP      | emit copy of last bit |
| 011  | FLIP     | emit complement of last bit |
| 100  | RPT2     | append 2 copies of last bit |
| 101  | RPT4     | append 4 copies |
| 110  | ALT      | append '01' |
| 111  | HALT     | stop |

Programs are 3-bit aligned bytestrings up to `L_MAX = 12` bits.  Programs
loop (PC wraps to 0 at end) so a 3-bit program "PUSH 0" emits an
arbitrarily long zero string with `Kt(0^n) = O(log n)`.

Bounded `Kt` for a target string is computed by full enumeration over
all programs of length at most `L_MAX` and minimising `|Pi| + ceil(log2 t)`
over Pi that match the target.  This is rich enough for:

- `Kt_bounded(0^n) = O(log n)` (compressible)
- `Kt_bounded((01)^n) = O(log n)` (compressible)
- `Kt_bounded(<random n-bit string>) = saturates at L_MAX + log T_MAX`
  (incompressible at our resolution)

This is sufficient to demonstrate Brandt's TRAVERSE in action and to
benchmark `pi_N` against random-baseline and against simple sequences.

## The MKtP query interface

`pi_as_MKtP_instance(N, k_offsets)` builds `(pi_N, k)` instances for k
ranging over `n + offset`.  This is what a hypothetical `Pi_Kt` solver
would be asked.  Note: this gives YES/NO answers about pi_N's Kt
*magnitude* but does NOT, even hypothetically, give us a TM that
produces a Kt-random string of growing length.  That gap is the
structural reason why Brandt's contradiction does not transfer
(see brandt_mktp_results.md).

## Edge IDs cited

- E5.3  (growing-dim MPOW in TC^0; the only open Chain-E sub-edge)
- E5.5  (Karchmer-Wigderson formula bound 2^{N/2-O(1)} — the only
          unconditional super-linear circuit-size lower bound on pi(x))
- E1.9  (pseudorandomness of pi(x) mod 2 under 33 measures)
- E7.10 (AKS modulus-twist orthogonality: Chain E "computationally cornered")
- T4    (Natural Proofs barrier: the constraint Brandt allegedly threads)
