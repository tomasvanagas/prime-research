/* Extended bounded-Kt VM simulator with arithmetic-primitive ablation.
 *
 * Identical to ../brandt_per_bit_arith_vm/sim.c except for an added
 * forbidden_ops_mask parameter on the batch evaluator: any program whose
 * nibble-decomposition contains a forbidden op is skipped (returned as
 * if it produced no output).  This lets us re-run the C3.a baseline with
 * any subset of {LOG2, LI, DIV_LOG, GEO_SUM} disabled to isolate which
 * primitives are causally responsible for the bit-J compression cut.
 *
 * Op codes (4-bit):  0 PUSH0, 1 PUSH1, 2 PUSH_N, 3 INC, 4 ADD, 5 SUB,
 *                    6 MUL,   7 SHR1,  8 LOG2,   9 LI, a DIV_LOG,
 *                    b GEO_SUM, c DUP, d EMIT_S, e EMIT0, f HALT.
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
    int sp = 0;
    uint64_t out_bits = 0;
    int out_count = 0;
    int steps = 0;
    int pc = 0;

    while (steps < T_MAX) {
        if (pc >= L_nyb) pc = 0;
        unsigned op = (prog_bits >> (4 * pc)) & 0xF;
        steps++;
        switch (op) {
            case 0xf: goto done;
            case 0x0: if (sp < STACK_CAP) S[sp++] = 0; break;
            case 0x1: if (sp < STACK_CAP) S[sp++] = 1; break;
            case 0x2: if (sp < STACK_CAP) S[sp++] = (int64_t)out_count; break;
            case 0x3:
                if (sp >= 1) {
                    int64_t v = S[sp - 1] + 1;
                    if (v > INT_CAP) v = INT_CAP;
                    S[sp - 1] = v;
                }
                break;
            case 0x4: {
                if (sp >= 2) {
                    int64_t a = S[--sp];
                    int64_t b = S[--sp];
                    int64_t v = a + b;
                    if (v > INT_CAP) v = INT_CAP;
                    S[sp++] = v;
                }
                break;
            }
            case 0x5: {
                if (sp >= 2) {
                    int64_t a = S[--sp];
                    int64_t b = S[--sp];
                    int64_t v = (b > a) ? (b - a) : 0;
                    S[sp++] = v;
                }
                break;
            }
            case 0x6: {
                if (sp >= 2) {
                    int64_t a = S[--sp];
                    int64_t b = S[--sp];
                    int64_t v = a * b;
                    if (v > INT_CAP || v < 0) v = INT_CAP;
                    S[sp++] = v;
                }
                break;
            }
            case 0x7:
                if (sp >= 1) S[sp - 1] >>= 1;
                break;
            case 0x8: {
                if (sp >= 1) {
                    int64_t v = S[sp - 1];
                    if (v < 1) S[sp - 1] = 0;
                    else {
                        int bl = 0;
                        int64_t t = v;
                        while (t) { bl++; t >>= 1; }
                        S[sp - 1] = bl - 1;
                    }
                }
                break;
            }
            case 0x9: {
                if (sp >= 1) {
                    int64_t a = S[--sp];
                    if (a < 2) S[sp++] = 0;
                    else {
                        double la = log((double)a);
                        int64_t v = (int64_t)((double)a / la);
                        if (v > INT_CAP) v = INT_CAP;
                        S[sp++] = v;
                    }
                }
                break;
            }
            case 0xa: {
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
            case 0xb: {
                if (sp >= 1) {
                    int64_t a = S[--sp];
                    if (a == 0) S[sp++] = 1;
                    else if (a == 1) S[sp++] = T_MAX;
                    else {
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
            case 0xc:
                if (sp >= 1 && sp < STACK_CAP) {
                    S[sp] = S[sp - 1];
                    sp++;
                }
                break;
            case 0xd: {
                int b = 0;
                if (sp >= 1) b = (int)(S[--sp] & 1);
                if (b) out_bits |= ((uint64_t)1 << out_count);
                out_count++;
                break;
            }
            case 0xe:
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

int emit_ops_in_cycle(uint64_t prog_bits, int L_nyb)
{
    int count = 0;
    for (int i = 0; i < L_nyb; i++) {
        unsigned op = (prog_bits >> (4 * i)) & 0xF;
        if (op == 0xf) return count;
        if (op == 0xd || op == 0xe) count++;
    }
    return count;
}

/* Returns true iff prog_bits contains any op whose bit is set in
 * forbidden_ops_mask.  Stops scanning at HALT. */
static inline int has_forbidden(uint64_t prog_bits, int L_nyb,
                                unsigned forbidden_ops_mask)
{
    if (forbidden_ops_mask == 0) return 0;
    for (int i = 0; i < L_nyb; i++) {
        unsigned op = (prog_bits >> (4 * i)) & 0xF;
        if (op == 0xf) return 0; /* HALT: rest never executes */
        if (forbidden_ops_mask & (1u << op)) return 1;
    }
    return 0;
}

/* Same as the baseline batch_enumerate_at_length but skips programs
 * that contain any op listed in forbidden_ops_mask (a 16-bit mask). */
int64_t batch_enumerate_at_length_filtered(
    int L_nyb,
    uint64_t *targets, int *target_lens, int n_targets,
    int *bests, uint64_t *best_progs, int *best_steps,
    int max_target_len,
    unsigned forbidden_ops_mask)
{
    int64_t n_progs = (int64_t)1 << (4 * L_nyb);
    int L_bits = 4 * L_nyb;
    int64_t simulated = 0;

    for (int64_t bits = 0; bits < n_progs; bits++) {
        if (emit_ops_in_cycle((uint64_t)bits, L_nyb) == 0) continue;
        if (has_forbidden((uint64_t)bits, L_nyb, forbidden_ops_mask)) continue;
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
