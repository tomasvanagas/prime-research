# wtrick_control.py — Results

Auxiliary script for the §C4 Anderson localisation attack. See the
primary writeup `anderson_localisation_chi_p_results.md` (Result 3) for
the full discussion.

**Summary:** the W-trick cascade. With W ∈ {6, 30, 210, 2310}, the
control potential is a random subset of `{n in [1, N] : gcd(n, W) = 1}`
of size matching chi_P's coprime-to-W count, plus fixed delta-functions
at the small primes p | W. Observed residual deviations
gamma_prime - gamma_CW (max over a 31-energy grid):

| W    | residual max |z|  | scale of W   |
|------|-------------------|--------------|
| 6    | 11.9 sigma        | mod 2, 3     |
| 30   | 6.3 sigma         | mod 2, 3, 5  |
| 210  | 6.1 sigma         | mod 2..7     |
| 2310 | 4.0 sigma         | mod 2..11    |

Cascade ends at borderline noise (~ 4 sigma on 31 energies, Bonferroni
threshold ~ 3.16). Output JSON: `wtrick_N{N}_s{seeds}_e{energies}.json`.

Confirms the spectral analogue of S85's Gowers-norm W-trick story.
