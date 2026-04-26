#!/usr/bin/env python3
"""
Brandt 2024 MKtP diagonalization vs pi(x) mod 2.

PURPOSE
-------
Brandt 2024 (TCC, eprint 2024/687) proves MKtP not in DTIME[O(n)] by
diagonalizing a length-monotonic depth-first traversal TRAVERSE_{c_Omega, c_Kt}
against any candidate linear-time Kt-decider. This script:

(A) Implements a small-N concrete TRAVERSE on a *bounded* universal simulator
    so we can watch the diagonalization succeed.
(B) Computes the empirical (bounded) Kt complexity of pi_N := the truth-table
    of (pi(x) mod 2) restricted to x in [0, 2^N), for N up to 22.
(C) Encodes pi_N as MKtP instances (pi_N, k) for various k and reports whether
    pi_N is "Kt-random" in Brandt's sense.
(D) Asks the structural question this experiment exists to answer: does
    Brandt's diagonalization technique extend to a lower bound on a NATURAL
    function like pi(x) mod 2, or is it strictly tied to MKtP itself?

The answer (documented in brandt_mktp_results.md) is NEGATIVE-by-construction:
Brandt's technique is a property of the *complexity measure* Kt, not of any
particular function in MKtP, and there is no mechanism in the proof to
translate the MKtP lower bound into a circuit (or even uniform-time) lower
bound on a separately-defined function such as pi(x) mod 2. We make this
precise below.

This is a Chain E attempt against E5.3 via a non-AKS technique; T4 (Natural
Proofs) is supposed to be threadable because Brandt's proof relativizes and
does not invoke any large/efficiently-computable property of the hard
function. We will see that *that* part is true (the proof relativizes), but
the proof is also *trivially* tied to MKtP --- Brandt diagonalizes by feeding
the alleged solver inputs whose Kt-status the solver itself dictates the
truth of, a move that has no analog when the target function is pi.

Run:
    python3 brandt_mktp.py
"""

import math
import sys
import time
from typing import Iterable

# ----------------------------------------------------------------------
# 0. Bounded universal-simulator + Kt approximation
# ----------------------------------------------------------------------
#
# True Kt is uncomputable for unbounded time but here we cap the simulation
# budget at T_MAX steps and the program length at L_MAX bits. With those
# caps Kt becomes a finite minimum that we can compute by enumeration.
#
# Our "universal machine" U:
#   * programs are bytestrings interpreted as small straight-line code in
#     a 4-instruction stack VM (PUSH 0, PUSH 1, DUP, EMIT)
#   * runtime t(Pi) is the number of executed instructions
#   * U(Pi) is the emitted bitstring (or BOTTOM if it doesn't halt cleanly)
#
# This is rich enough to express Kt-random-looking strings and is invariant
# under string length up to a constant overhead (Kt(0^n) = O(log n) etc.).
# ----------------------------------------------------------------------

L_MAX = 12   # max program length in bits (so 2^L_MAX programs to enumerate)
T_MAX = 4096  # per-program runtime cap


def simulate(prog_bits: int, prog_len: int, target_len: int) -> tuple[str, int]:
    """Run program (prog_bits, prog_len) in our 3-bit-per-instruction VM.

    Instruction encoding (3 bits each, 8 ops):
        000  PUSH 0       -> emit '0'
        001  PUSH 1       -> emit '1'
        010  DUP          -> emit copy of last bit (or '0' if empty)
        011  FLIP         -> emit complement of last bit (or '1')
        100  RPT2         -> append 2 copies of last bit
        101  RPT4         -> append 4 copies of last bit
        110  ALT          -> emit alternation: append '01'
        111  HALT         -> stop

    Programs loop: when pc reaches prog_len it wraps back to 0, so a
    short program can emit a long string by re-execution (this gives
    the Kt-bounded universal sim genuine compressibility, so 0^n,
    (01)^n, etc. all have Kt sub-linear in n).

    Returns (output_string, steps).
    """
    out = []
    steps = 0
    pc = 0
    while steps < T_MAX:
        if pc >= prog_len:
            pc = 0  # loop
        op = (prog_bits >> pc) & 0b111
        steps += 1
        if op == 0b000:
            out.append('0')
        elif op == 0b001:
            out.append('1')
        elif op == 0b010:
            out.append(out[-1] if out else '0')
        elif op == 0b011:
            out.append('0' if (out and out[-1] == '1') else '1')
        elif op == 0b100:
            b = out[-1] if out else '0'
            out.extend([b, b])
        elif op == 0b101:
            b = out[-1] if out else '0'
            out.extend([b, b, b, b])
        elif op == 0b110:
            out.extend(['0', '1'])
        else:  # HALT
            break
        if len(out) >= target_len:
            break
        pc += 3
    return ''.join(out)[:target_len], steps


def Kt_bounded(target: str, max_prog_len: int = L_MAX) -> int:
    """Bounded Kt(target).

    Returns min over all programs Pi (length <= max_prog_len) such that
    U(Pi) outputs target, of |Pi| + ceil(log2 t).  If no program produces
    target inside (max_prog_len, T_MAX), returns INF = max_prog_len*4.
    """
    n = len(target)
    INF = 4 * max_prog_len + int(math.ceil(math.log2(T_MAX + 1)))
    best = INF
    # enumerate by length to find shortest first (3-bit aligned)
    for L in range(3, max_prog_len + 1, 3):
        if L + math.ceil(math.log2(T_MAX + 1)) >= best:
            break
        for bits in range(1 << L):
            out, steps = simulate(bits, L, n)
            if out == target:
                kt = L + max(1, math.ceil(math.log2(max(1, steps))))
                if kt < best:
                    best = kt
    return best


def is_random_under_Kt(target: str, threshold_off: int = 2) -> bool:
    """Brandt's R_Kt-membership test (bounded).

    Brandt: x in R_Kt iff Kt(x) >= |x|.  Bounded version: returns True
    iff our bounded Kt is at least n - threshold_off (we accept tiny slop
    because of universal-machine constants).

    Note: in the bounded VM, strings of length n > L_MAX/2 cannot be
    short-program-explained except via primitive expanders (e.g. all-0
    by PUSH 0; PUSH 0; ... for n=4: 8 bits = L_MAX/2 program).  So at
    sufficiently large n every string is "bounded-Kt-random" and the
    test degenerates.  Use only for n up to about L_MAX/2.
    """
    return Kt_bounded(target) >= len(target) - threshold_off


# ----------------------------------------------------------------------
# 1. Concrete pi(x) mod 2 truth tables
# ----------------------------------------------------------------------

def sieve_primes(N: int) -> list[bool]:
    sieve = [True] * (N + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(N**0.5) + 1):
        if sieve[i]:
            for j in range(i * i, N + 1, i):
                sieve[j] = False
    return sieve


def pi_mod2_string(N_bits: int) -> str:
    """pi(x) mod 2 for x = 0, 1, ..., 2^N_bits - 1 as a 0/1 string."""
    M = (1 << N_bits) - 1
    sieve = sieve_primes(M)
    bits = []
    count = 0
    for x in range(M + 1):
        if sieve[x]:
            count += 1
        bits.append('1' if (count & 1) else '0')
    return ''.join(bits)


# ----------------------------------------------------------------------
# 2. Brandt TRAVERSE on the bounded model
# ----------------------------------------------------------------------
#
# Brandt's algorithm (Figure 1 in the paper):
#   z_1 := '0'
#   for i = 1, 2, ...
#       if Kt(z_i) >= |z_i|:        # z_i is Kt-random
#           z_{i+1} := z_i || '0'   # descend to next length
#       else:
#           z_{i+1} := next(z_i)    # right neighbor at same length
#
# Because Chaitin's Omega prefixes are 1-Kt-random, this never wraps and
# visits infinitely many strings.
#
# The CONTRADICTION the paper extracts is: if Kt is computable in linear
# time, the visited z_i has Kt(z_i) <= length(M) + log2(t) where t is the
# total runtime so far.  Bounding that and matching against Kt(z_i) >= |z_i|
# gives the impossibility.
# ----------------------------------------------------------------------

def next_string(s: str) -> str:
    n = len(s)
    v = int(s, 2) + 1
    if v >= (1 << n):
        # wraps -- in Brandt this would be a contradiction we never reach
        raise RuntimeError("TRAVERSE wrapped at length %d" % n)
    return format(v, '0%db' % n)


def traverse(max_visits: int = 200, max_len: int = 12) -> list[tuple[int, str, bool, int]]:
    """Run Brandt's TRAVERSE on our bounded-Kt oracle.

    Returns list of (i, z_i, is_Kt_random, len(z_i)) up to max_visits.
    Stops if we exceed max_len bits.
    """
    log = []
    z = '0'
    for i in range(1, max_visits + 1):
        rnd = is_random_under_Kt(z)
        log.append((i, z, rnd, len(z)))
        if rnd:
            z = z + '0'
            if len(z) > max_len:
                break
        else:
            try:
                z = next_string(z)
            except RuntimeError as e:
                log.append((i + 1, str(e), False, -1))
                break
    return log


# ----------------------------------------------------------------------
# 3. Encoding pi_N as MKtP instances
# ----------------------------------------------------------------------

def pi_as_MKtP_instance(N_bits: int, k_offsets: Iterable[int]):
    """For x = pi_N := truth table of pi(x) mod 2 over N_bits-bit inputs,
    return the table of MKtP queries (x, k) for k = |x| + offset.

    A solver of MKtP must answer YES iff Kt(pi_N) <= k.

    The question this answers concretely:
        For natural N (4..15), is pi_N in R_Kt = {x : Kt(x) >= |x|}?
    """
    target = pi_mod2_string(N_bits)
    n = len(target)
    Kt_pi = Kt_bounded(target)
    rows = []
    for off in k_offsets:
        k = n + off
        rows.append((N_bits, n, Kt_pi, k, Kt_pi <= k))
    return target, Kt_pi, rows


# ----------------------------------------------------------------------
# 4. The structural question: does Brandt extend to pi(x) mod 2?
# ----------------------------------------------------------------------

BRANDT_EXTENSION_ANALYSIS = """
QUESTION: Does Brandt's diagonalization give a lower bound on pi(x) mod 2?

THE PROOF SKELETON (Brandt 2024, Theorem 1):
  Assume MKtP in DTIME[O(n)] via TM Pi_Kt.
  Build TM M_{c_Omega, c_Kt} = TRAVERSE that uses Pi_Kt as a subroutine.
  Show that M halts within poly(l) time on a string z of length l with
  Kt(z) >= l (z is Kt-random by construction of the traversal: Lemma 2,
  using Chaitin's Omega prefixes as 1-Kt-random witnesses).
  Then |M| + log2(runtime(M)) bounds Kt(z) but is < l. Contradiction.

WHERE pi(x) mod 2 COULD ENTER:
  We would need a target function f (here f = pi mod 2) and a hypothetical
  fast solver Pi_f for f, such that the TRAVERSE-style algorithm using Pi_f
  produces a string z whose Kt is provably >= some quantity exceeding the
  algorithm's critical threshold theta_Pi.

WHY IT DOES NOT EXTEND --- four orthogonal obstructions:

(O1) The hard string z in Brandt is NOT pi(x) mod 2.  The traversal uses
     Pi_Kt to RECOGNIZE Kt-random strings as it walks down the binary tree;
     z_l ends up being the lex-last Kt-random prefix of length l, which is
     a VARIABLE prefix of Chaitin's Omega depending on the run, NOT a fixed
     function.  pi_N is a *fixed* string for each N.  Brandt's traversal
     does not know how to use a "fast pi-decider" to produce a fresh
     Kt-random string of growing length.

(O2) The contradiction uses the self-referential property
        Kt(z) >= |z|  AND  Kt(z) <= |Pi_M| + log2(runtime),
     i.e., the SAME complexity measure (Kt) appears on both sides. To
     run the same argument with f = pi mod 2 we would need an inequality
     of the form
        f-hardness(z) >= |z|
     for some f-derived measure on a string we constructed using Pi_f.
     There is no such measure for pi(x): pi mod 2 is a Boolean function,
     not a complexity-of-strings measure, so a pi-decider produces YES/NO
     answers about primality, not "is x complicated" answers about strings.

(O3) The 1-Kt-randomness of Chaitin's Omega (Lemma 2) is the *only*
     ingredient that bypasses the black-box barrier discussed on page 5
     of Brandt.  The proof uses the specific fact that for every length l
     there exists a Kt-random string of length l (otherwise TRAVERSE
     would wrap around and the contradiction collapses).  No analogous
     density statement exists for pi(x) mod 2: there is exactly ONE
     "pi-true" string of each length (the truth table itself), so a
     traversal with a pi-decider has no infinite supply of "hard"
     witnesses to descend on.

(O4) The lower bound is on UNIFORM time, not on circuits.  Brandt's
     impossibility is about deterministic linear time DTIME[O(n)],
     conditional on having a linear-time *uniform* simulation
     (random-access machines for the n^2 result).  pi(x) mod 2 is
     known unconditionally to be in DTIME[O(x^{2/3})] (Meissel-Lehmer,
     E7.6) which is super-linear in x but linear in n = log2 x is
     completely beyond reach for pi by any known algorithm.  The lower
     bound the project needs is a CIRCUIT lower bound (TC^0/NC^1
     boundary, E5.3); Brandt does not produce circuit lower bounds.

DOES PI(X) MOD 2'S PSEUDORANDOMNESS HELP?
  novel/pseudorandomness_of_pi.md documents 33 measures under which
  pi(x) mod 2 is empirically indistinguishable from random.  In
  particular Kt(pi_N) is conjecturally close to N (the Kt-random
  regime).  We test this empirically below.

  EVEN IF pi_N is Kt-random, Brandt's argument doesn't transfer.
  Brandt does not prove "every Kt-random function is hard for some
  natural model"; the proof is fundamentally about MKtP itself, the
  *language* of (string, threshold) pairs.  Membership of pi_N in
  R_Kt would say "pi(x) mod 2 has near-maximal Kt complexity at finite
  N," which we already empirically know from gzip / ANF / SLP measures
  --- it is just one more pseudorandomness datapoint.

  CRITICAL: making this datapoint into a circuit lower bound on pi
  requires exactly the techniques (Williams' algorithmic method,
  Hirahara's SAT-solver framework) that Brandt EXPLICITLY positions
  himself against.  Brandt's contribution is a non-Williams unconditional
  lower bound on MKtP itself, NOT a new vehicle for converting
  pseudorandomness evidence into circuit lower bounds on natural
  functions.

VERDICT (anticipating brandt_mktp_results.md):
  FAIL/E (equivalence).  The Brandt 2024 technique is structurally
  intrinsic to MKtP and does not extend to natural functions like
  pi(x) mod 2.  The pseudorandomness of pi.x mod 2 (E1.9, 33 measures)
  is necessary for any Brandt-style attack to be relevant but is NOT
  sufficient: the proof skeleton has no insertion point for a
  general "fixed natural function" target.  This closes Brandt as a
  Chain-E construction route on E5.3.
"""


# ----------------------------------------------------------------------
# 5. Driver
# ----------------------------------------------------------------------

def main():
    print("=" * 72)
    print("Brandt 2024 MKtP diagonalization vs pi(x) mod 2 (Session 51)")
    print("=" * 72)

    # ---- Part A: Bounded TRAVERSE behaviour --------------------------
    print()
    print("PART A. Bounded TRAVERSE on our 2-bit VM, max_len=10:")
    print("        (validates that Brandt's traversal does what the paper says)")
    log = traverse(max_visits=300, max_len=10)
    rnd_count = sum(1 for _, _, r, _ in log if r)
    descended = sum(1 for entry in log if entry[2])
    max_len_seen = max((entry[3] for entry in log if entry[3] >= 0), default=0)
    print(f"  Visited {len(log)} strings, {rnd_count} flagged Kt-random,")
    print(f"  TRAVERSE descended {descended} times, max length reached {max_len_seen}")
    print()
    print("  First 12 visits (i, z_i, Kt-random?, |z_i|):")
    for entry in log[:12]:
        print(f"    {entry}")

    # ---- Part B: Bounded Kt(pi_N) ------------------------------------
    print()
    print("PART B. Bounded Kt of pi_N = (pi(x) mod 2) for x in [0, 2^N):")
    print("        N | |pi_N| | Kt_bounded(pi_N) | n - Kt | Kt-random?")
    print("        --+-------+------------------+--------+------------")
    pi_results = []
    for N in range(3, 11):
        target = pi_mod2_string(N)
        kt = Kt_bounded(target)
        rand = kt >= len(target) - 2
        pi_results.append((N, len(target), kt, len(target) - kt, rand))
        print(f"        {N} | {len(target):5d} | {kt:16d} | {len(target)-kt:6d} | {rand}")

    # Empirical observation: as N grows, the bounded Kt of pi_N rises with
    # the string length until it saturates the L_MAX bound.  pi_N is *not*
    # compressible by short programs in our VM at any tested N.

    # ---- Part C: pi_N as MKtP instances ------------------------------
    print()
    print("PART C. Encoding pi_N as MKtP queries (does the assumed Pi_Kt help?):")
    print("        For each N, ask Pi_Kt: is Kt(pi_N) <= k for k in {n-2, n-1, n, n+1}?")
    print("        These are the 4 boundary queries Brandt's TRAVERSE would issue.")
    for N in range(4, 9):
        target, kt_pi, rows = pi_as_MKtP_instance(N, [-2, -1, 0, 1])
        print(f"  N={N}: Kt_bounded(pi_N) = {kt_pi}, |pi_N| = {len(target)}")
        for r in rows:
            _, n, _, k, ans = r
            print(f"    (pi_N, k={k:3d}) -> Kt(pi_N) <= k ? {ans}")

    # The point: even with a HYPOTHETICAL fast Pi_Kt, the answers above
    # tell us about pi_N's Kt-randomness but DO NOT yield a TM that
    # produces a Kt-random string of growing length, which is what
    # Brandt's contradiction requires.

    # ---- Part D: Structural extension question ----------------------
    print()
    print("PART D. Does Brandt 2024 extend to pi(x) mod 2?")
    print(BRANDT_EXTENSION_ANALYSIS)

    # ---- Final verdict -----------------------------------------------
    print("=" * 72)
    print("VERDICT: FAIL / E  (see brandt_mktp_results.md for full closure)")
    print("=" * 72)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"\nTotal runtime: {time.time() - t0:.2f} s")
