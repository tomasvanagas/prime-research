/* Extended bounded-Kt VM simulator (C implementation for ctypes).
 *
 * Op codes (4-bit):
 *   0  PUSH0       push 0 onto integer stack
 *   1  PUSH1       push 1
 *   2  PUSH_N      push current emit-count (out_count)
 *   3  INC         TOS += 1 (saturating at INT_CAP)
 *   4  ADD         pop a, b; push min(a+b, INT_CAP)
 *   5  SUB         pop a, b; push max(b-a, 0)
 *   6  MUL         pop a, b; push min(a*b, INT_CAP)
 *   7  SHR1        TOS >>= 1
 *   8  LOG2        TOS <- floor(log2(TOS))   (log2(0) := 0)
 *   9  LI          pop a; push floor(a / log a) for a >= 2 else 0
 *   a  DIV_LOG     pop a, b; push floor(b / log max(a, 2))
 *   b  GEO_SUM     pop a; push 1 + a + a^2 + ... while sum <= T_MAX
 *   c  DUP         duplicate TOS
 *   d  EMIT_S      pop a; emit (a & 1) bit
 *   e  EMIT0       emit 0 bit
 *   f  HALT        stop
 *
 * Programs of length L (in nibbles) loop: pc wraps to 0 when reaching L.
 *
 * Returns the count of bits emitted (<= target_len).  Sets *out_bits to
 * the LSB-first packed bit pattern, *steps_out to the number of executed
 * instructions when termination occurred.
 */

#include <stdint.h>
#include <math.h>
#include <stdbool.h>

#define T_MAX  4096
#define INT_CAP 1000000000LL
#define STACK_CAP 32

int simulate(
    uint64_t prog_bits, int L_nyb, int target_len,
    uint64_t *out_bits_out, int *steps_out)
{
    int64_t S[STACK_CAP];
    int sp = 0;            /* stack pointer (number of items) */
    uint64_t out_bits = 0;
    int out_count = 0;
    int steps = 0;
    int pc = 0;

    while (steps < T_MAX) {
        if (pc >= L_nyb) pc = 0;
        unsigned op = (prog_bits >> (4 * pc)) & 0xF;
        steps++;
        switch (op) {
            case 0xf:  /* HALT */
                goto done;
            case 0x0:  /* PUSH0 */
                if (sp < STACK_CAP) S[sp++] = 0;
                break;
            case 0x1:  /* PUSH1 */
                if (sp < STACK_CAP) S[sp++] = 1;
                break;
            case 0x2:  /* PUSH_N */
                if (sp < STACK_CAP) S[sp++] = (int64_t)out_count;
                break;
            case 0x3:  /* INC */
                if (sp >= 1) {
                    int64_t v = S[sp - 1] + 1;
                    if (v > INT_CAP) v = INT_CAP;
                    S[sp - 1] = v;
                }
                break;
            case 0x4: { /* ADD */
                if (sp >= 2) {
                    int64_t a = S[--sp];
                    int64_t b = S[--sp];
                    int64_t v = a + b;
                    if (v > INT_CAP) v = INT_CAP;
                    S[sp++] = v;
                }
                break;
            }
            case 0x5: { /* SUB */
                if (sp >= 2) {
                    int64_t a = S[--sp];
                    int64_t b = S[--sp];
                    int64_t v = (b > a) ? (b - a) : 0;
                    S[sp++] = v;
                }
                break;
            }
            case 0x6: { /* MUL */
                if (sp >= 2) {
                    int64_t a = S[--sp];
                    int64_t b = S[--sp];
                    int64_t v = a * b;
                    if (v > INT_CAP || v < 0) v = INT_CAP;
                    S[sp++] = v;
                }
                break;
            }
            case 0x7:  /* SHR1 */
                if (sp >= 1) {
                    S[sp - 1] >>= 1;
                }
                break;
            case 0x8: { /* LOG2 */
                if (sp >= 1) {
                    int64_t v = S[sp - 1];
                    if (v < 1) {
                        S[sp - 1] = 0;
                    } else {
                        /* Find bit_length - 1 */
                        int bl = 0;
                        int64_t t = v;
                        while (t) { bl++; t >>= 1; }
                        S[sp - 1] = bl - 1;
                    }
                }
                break;
            }
            case 0x9: { /* LI */
                if (sp >= 1) {
                    int64_t a = S[--sp];
                    if (a < 2) {
                        S[sp++] = 0;
                    } else {
                        double la = log((double)a);
                        int64_t v = (int64_t)((double)a / la);
                        if (v > INT_CAP) v = INT_CAP;
                        S[sp++] = v;
                    }
                }
                break;
            }
            case 0xa: { /* DIV_LOG */
                if (sp >= 2) {
                    int64_t a = S[--sp];
                    int64_t b = S[--sp];
                    int64_t aa = (a >= 2) ? a : 2;
                    double la = log((double)aa);
                    int64_t v = (int64_t)((double)b / la);
                    if (v > INT_CAP) v = INT_CAP;
                    S[sp++] = v;
                }
                break;
            }
            case 0xb: { /* GEO_SUM */
                if (sp >= 1) {
                    int64_t a = S[--sp];
                    if (a == 0) {
                        S[sp++] = 1;
                    } else if (a == 1) {
                        S[sp++] = T_MAX;
                    } else {
                        int64_t s = 1, p = 1;
                        while (1) {
                            p *= a;
                            if (p < 0 || s + p > T_MAX) break;
                            s += p;
                            if (s >= INT_CAP) { s = INT_CAP; break; }
                        }
                        S[sp++] = s;
                    }
                }
                break;
            }
            case 0xc:  /* DUP */
                if (sp >= 1 && sp < STACK_CAP) {
                    S[sp] = S[sp - 1];
                    sp++;
                }
                break;
            case 0xd: { /* EMIT_S */
                int b = 0;
                if (sp >= 1) {
                    b = (int)(S[--sp] & 1);
                }
                if (b) out_bits |= ((uint64_t)1 << out_count);
                out_count++;
                break;
            }
            case 0xe:  /* EMIT0 */
                out_count++;
                break;
        }
        if (out_count >= target_len) goto done;
        pc++;
    }
done:
    *out_bits_out = out_bits;
    *steps_out = steps;
    return out_count;
}

/* Pre-check: returns the number of EMIT ops in one pc cycle, stopping at
 * HALT.  Used to skip programs that produce no output. */
int emit_ops_in_cycle(uint64_t prog_bits, int L_nyb)
{
    int count = 0;
    for (int i = 0; i < L_nyb; i++) {
        unsigned op = (prog_bits >> (4 * i)) & 0xF;
        if (op == 0xf) return count;             /* HALT */
        if (op == 0xd || op == 0xe) count++;     /* EMIT */
    }
    return count;
}

/* Batch evaluator: enumerate every program of length L_nyb (in nibbles),
 * skip programs with 0 emits in cycle, and for each candidate that emits
 * MAX_TARGET_LEN bits (or hits the requested target_len), check against
 * each of n_targets target patterns.  Updates *bests in-place.
 *
 * targets[i]   : target bit pattern (LSB-first packed, length target_lens[i])
 * target_lens  : individual lengths
 * n_targets    : number of targets
 * bests        : current best Kt for each target (in/out)
 * best_progs   : program bits achieving bests[i]
 * best_steps   : steps used in that simulation
 *
 * For each program: simulate to max(target_lens).  For each target, check
 * if its prefix matches.  If yes and current Kt = L_bits + ceil(log2(steps))
 * is below bests[i], update.
 *
 * Returns the number of programs actually simulated (skipped programs not
 * counted).
 */
int64_t batch_enumerate_at_length(
    int L_nyb,
    uint64_t *targets, int *target_lens, int n_targets,
    int *bests, uint64_t *best_progs, int *best_steps,
    int max_target_len)
{
    int64_t n_progs = (int64_t)1 << (4 * L_nyb);
    int L_bits = 4 * L_nyb;
    int64_t simulated = 0;

    for (int64_t bits = 0; bits < n_progs; bits++) {
        if (emit_ops_in_cycle((uint64_t)bits, L_nyb) == 0) continue;
        uint64_t out_bits;
        int steps;
        int oc = simulate((uint64_t)bits, L_nyb, max_target_len,
                          &out_bits, &steps);
        simulated++;
        for (int i = 0; i < n_targets; i++) {
            int tl = target_lens[i];
            if (oc < tl) continue;
            uint64_t mask = (tl >= 64) ? ~(uint64_t)0 : (((uint64_t)1 << tl) - 1);
            if ((out_bits & mask) == targets[i]) {
                /* compute kt */
                int log_steps = 1;
                int s = steps;
                if (s > 1) {
                    log_steps = 0;
                    int s2 = s - 1;
                    while (s2) { log_steps++; s2 >>= 1; }
                }
                int kt = L_bits + log_steps;
                if (kt < bests[i]) {
                    bests[i] = kt;
                    best_progs[i] = (uint64_t)bits;
                    best_steps[i] = steps;
                }
            }
        }
    }
    return simulated;
}
