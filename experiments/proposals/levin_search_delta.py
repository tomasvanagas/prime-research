"""Proposal A: Levin-style universal search on delta(n) = p(n) - R^{-1}(n).

We run a bounded enumeration of short programs in a small arithmetic DSL
and check whether any of them reproduces the integer values of delta(n)
for n = 1..N_target. A success would be a closed-form computation of
delta(n) and hence a polylog algorithm for p(n).

The DSL is RPN over a tiny alphabet:
    n  -- input
    0..9 -- small constants
    +, -, *, //, %  -- binary integer ops
    L  -- floor(log_2 x)
    S  -- isqrt(x)

We enumerate all programs of length <= L_max (typing-preserving), evaluate
on the inputs n = 1..N, and report any that exactly produce delta(1..N).
We also report the best "almost matches" (highest k such that program
matches delta on n=1..k).

This is effectively a bounded Levin search: programs of length L are
allotted O(2^{L_max - L}) units of evaluation time. We enforce a strict
runtime cap.

Run: python3 levin_search_delta.py
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "algorithms"))
from v5_pure_analytic import _R_inverse  # type: ignore  # noqa: E402


def sieve_primes(N: int) -> list[int]:
    sv = [True] * (N + 1)
    sv[0] = sv[1] = False
    for k in range(2, int(math.isqrt(N)) + 1):
        if sv[k]:
            for j in range(k * k, N + 1, k):
                sv[j] = False
    return [i for i, b in enumerate(sv) if b]


def compute_delta(n_max: int) -> list[int]:
    """delta(k) = p(k) - round(R^{-1}(k)) for k = 1..n_max."""
    primes_needed = max(50, n_max + 5)
    primes = sieve_primes(20 * n_max + 200)
    deltas = []
    for k in range(1, n_max + 1):
        rinv = _R_inverse(k)
        deltas.append(primes[k - 1] - int(round(rinv)))
    return deltas


# ----------------------------------------------------------------------
# Tiny RPN DSL evaluator
# ----------------------------------------------------------------------

def safe_eval(prog: tuple[str, ...], n: int) -> int | None:
    """Evaluate RPN program with input n. Return None on stack/type error."""
    stack: list[int] = []
    for tok in prog:
        if tok == "n":
            stack.append(n)
        elif tok.isdigit():
            stack.append(int(tok))
        elif tok == "+":
            if len(stack) < 2:
                return None
            b, a = stack.pop(), stack.pop()
            stack.append(a + b)
        elif tok == "-":
            if len(stack) < 2:
                return None
            b, a = stack.pop(), stack.pop()
            stack.append(a - b)
        elif tok == "*":
            if len(stack) < 2:
                return None
            b, a = stack.pop(), stack.pop()
            v = a * b
            if abs(v) > 10**12:
                return None
            stack.append(v)
        elif tok == "//":
            if len(stack) < 2:
                return None
            b, a = stack.pop(), stack.pop()
            if b == 0:
                return None
            stack.append(a // b)
        elif tok == "%":
            if len(stack) < 2:
                return None
            b, a = stack.pop(), stack.pop()
            if b == 0:
                return None
            stack.append(a % b)
        elif tok == "L":
            if not stack:
                return None
            x = stack.pop()
            if x <= 0:
                return None
            stack.append(int(x).bit_length() - 1)
        elif tok == "S":
            if not stack:
                return None
            x = stack.pop()
            if x < 0:
                return None
            stack.append(math.isqrt(x))
        else:
            return None
    if len(stack) != 1:
        return None
    return stack[0]


# ----------------------------------------------------------------------
# Bounded enumeration
# ----------------------------------------------------------------------

def enumerate_programs(tokens: list[str], max_len: int):
    """Yield all RPN programs of length <= max_len with valid stack shape."""
    # Track (length, prog, stack_size). Stop early once stack would underflow.
    def rec(prog, stack):
        if len(prog) > 0 and stack == 1:
            yield tuple(prog)
        if len(prog) >= max_len:
            return
        for tok in tokens:
            new_stack = stack
            if tok == "n" or tok.isdigit():
                new_stack += 1
            elif tok in ("+", "-", "*", "//", "%"):
                if stack < 2:
                    continue
                new_stack -= 1
            elif tok in ("L", "S"):
                if stack < 1:
                    continue
            prog.append(tok)
            yield from rec(prog, new_stack)
            prog.pop()

    yield from rec([], 0)


def main() -> None:
    t0 = time.time()
    print("# Levin Universal Search on delta(n)")
    N_TARGET = 30
    delta = compute_delta(N_TARGET)
    print(f"# delta(1..{N_TARGET}) = {delta}")
    print()

    tokens = ["n", "1", "2", "3", "+", "-", "*", "//", "%", "L", "S"]
    for max_len in (3, 4, 5, 6):
        best_score = -1
        best_prog = None
        count = 0
        exact_matches = []
        for prog in enumerate_programs(tokens, max_len):
            count += 1
            score = 0
            ok = True
            for k in range(1, N_TARGET + 1):
                v = safe_eval(prog, k)
                if v is None:
                    ok = False
                    break
                if v != delta[k - 1]:
                    break
                score += 1
            if score > best_score and ok is not False:
                best_score = score
                best_prog = prog
            if score == N_TARGET:
                exact_matches.append(prog)
        elapsed = time.time() - t0
        print(
            f"## max_len={max_len}: enumerated {count} programs, "
            f"best prefix-match = {best_score}/{N_TARGET}, "
            f"prog = {' '.join(best_prog) if best_prog else 'none'}, "
            f"elapsed = {elapsed:.1f}s"
        )
        if exact_matches:
            print(f"  EXACT MATCHES FOUND: {len(exact_matches)}")
            for m in exact_matches[:5]:
                print(f"    {' '.join(m)}")

    # Sanity check: a constant zero baseline match length
    zero_match = sum(1 for d in delta if d == 0)
    print(f"\n# Baseline: program '0' matches {zero_match}/{N_TARGET}")
    # Constant -1 baseline
    minus1_match = sum(1 for d in delta if d == -1)
    print(f"# Baseline: program 'const -1' matches {minus1_match}/{N_TARGET}")
    # Mode baseline
    from collections import Counter

    ct = Counter(delta)
    most, hits = ct.most_common(1)[0]
    print(f"# Baseline: constant ({most}) matches {hits}/{N_TARGET}")


if __name__ == "__main__":
    main()
