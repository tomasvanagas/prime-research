"""
DARTS-style differentiable architecture search for TC^0 PRIMES circuits.

Attack vector: ATTACK_VECTORS.md §D8.
Cross-domain technique: Liu-Simonyan-Yang 2019 "DARTS: Differentiable
Architecture Search" (ICLR), Bender et al. 2018 "Understanding and
Simplifying One-Shot Architecture Search" (ICML).

Search space:
  Depth-3 circuit on N-bit inputs.
  Layer 1: G1 nodes; Layer 2: G2 nodes; Layer 3: 1 output.
  Each node is a softmax mixture over the gate library
    {AND, OR, XOR, MAJ_3, MAJ_5, MAJ_7, ID, NOT}
  with per-input soft selection logits beta in R^{prev_size}.

Continuous-relaxation gate definitions (sigmoid s_i := sigmoid(beta_i)):
  AND: prod_i (1 - s_i * (1 - v_i))
  OR : 1 - prod_i (1 - s_i * v_i)
  XOR: (1 - prod_i (1 - 2 s_i v_i)) / 2
  MAJ_k: sigmoid(tau * (sum_i s_i * v_i - k/2))   (k = 3, 5, 7)
  ID : s_first * v_first + (1 - s_first) * 0.5
  NOT: s_first * (1 - v_first) + (1 - s_first) * 0.5

Targets:
  - PRIMES: chi_P(n) = 1 if n prime else 0, for n in [0, 2^N)
  - Density-matched random Boolean function (control)

Evaluation:
  - Final BCE loss on training set
  - Discretised circuit accuracy
  - Generalisation to held-out n in [2^N, 2^N + 1000)
"""
import argparse
import json
import os
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from sympy import isprime

torch.set_num_threads(max(1, os.cpu_count() // 4))
torch.set_num_interop_threads(1)


def chi_P(n_array: np.ndarray) -> np.ndarray:
    """Vectorised primality indicator."""
    out = np.zeros_like(n_array, dtype=np.float32)
    for i, n in enumerate(n_array.tolist()):
        out[i] = 1.0 if (n >= 2 and isprime(int(n))) else 0.0
    return out


def n_to_bits(n_array: np.ndarray, N: int) -> np.ndarray:
    """Little-endian bit decomposition of a 1-D int array to shape (len, N)."""
    bits = ((n_array.reshape(-1, 1) >> np.arange(N).reshape(1, N)) & 1).astype(np.float32)
    return bits


GATE_LIBRARY = ["AND", "OR", "XOR", "MAJ3", "MAJ5", "MAJ7", "ID", "NOT"]
N_GATES = len(GATE_LIBRARY)


class DartsNode(nn.Module):
    """One differentiable node: softmax-mixture of gates over soft-selected inputs."""

    def __init__(self, prev_size: int, init_tau: float = 4.0):
        super().__init__()
        self.prev_size = prev_size
        # Architecture logits (softmax over gates).
        self.alpha = nn.Parameter(torch.zeros(N_GATES))
        # Input-selection logits (sigmoid for soft selection).
        self.beta = nn.Parameter(torch.zeros(prev_size))
        # Gate temperature (per-node, learned).
        self.log_tau = nn.Parameter(torch.tensor(np.log(init_tau), dtype=torch.float32))
        # First-input selection logit for ID/NOT gates.
        self.gamma = nn.Parameter(torch.zeros(prev_size))

    def forward(self, V: torch.Tensor) -> torch.Tensor:
        """V: (batch, prev_size) values in [0,1]. Returns (batch,) values in [0,1]."""
        s = torch.sigmoid(self.beta)            # (prev_size,)
        s_b = s.unsqueeze(0)                    # (1, prev_size)
        eps = 1e-6
        v_clamped = V.clamp(eps, 1 - eps)

        # AND: prod_i (1 - s_i (1 - v_i))
        and_terms = 1.0 - s_b * (1.0 - v_clamped)
        log_and = torch.log(and_terms.clamp(min=eps)).sum(dim=1)
        y_and = torch.exp(log_and)

        # OR: 1 - prod_i (1 - s_i v_i)
        or_terms = 1.0 - s_b * v_clamped
        log_or_neg = torch.log(or_terms.clamp(min=eps)).sum(dim=1)
        y_or = 1.0 - torch.exp(log_or_neg)

        # XOR: (1 - prod_i (1 - 2 s_i v_i))/2
        xor_terms = 1.0 - 2.0 * s_b * v_clamped       # in [-1, 1]
        # Use abs+sign trick to take log of negative-allowed product.
        sign = torch.sign(xor_terms).prod(dim=1)
        log_abs = torch.log(xor_terms.abs().clamp(min=eps)).sum(dim=1)
        prod_xor = sign * torch.exp(log_abs)
        y_xor = (1.0 - prod_xor) / 2.0

        # MAJ_k: sigmoid(tau (sum s_i v_i - k/2))
        tau = torch.exp(self.log_tau)
        weighted_sum = (s_b * v_clamped).sum(dim=1)         # (batch,)
        y_maj3 = torch.sigmoid(tau * (weighted_sum - 1.5))
        y_maj5 = torch.sigmoid(tau * (weighted_sum - 2.5))
        y_maj7 = torch.sigmoid(tau * (weighted_sum - 3.5))

        # ID / NOT — soft pick of one input by gamma.
        gamma_w = F.softmax(self.gamma, dim=0).unsqueeze(0)  # (1, prev_size)
        v_pick = (gamma_w * v_clamped).sum(dim=1)            # (batch,)
        y_id = v_pick
        y_not = 1.0 - v_pick

        outputs = torch.stack([y_and, y_or, y_xor, y_maj3, y_maj5, y_maj7, y_id, y_not], dim=1)
        # outputs: (batch, N_GATES)
        weights = F.softmax(self.alpha, dim=0).unsqueeze(0)  # (1, N_GATES)
        return (weights * outputs).sum(dim=1)


class DartsCircuit(nn.Module):
    def __init__(self, N: int, G1: int = 16, G2: int = 16, init_tau: float = 4.0):
        super().__init__()
        self.N = N
        self.G1 = G1
        self.G2 = G2
        self.layer1 = nn.ModuleList([DartsNode(prev_size=N, init_tau=init_tau) for _ in range(G1)])
        self.layer2 = nn.ModuleList([DartsNode(prev_size=G1, init_tau=init_tau) for _ in range(G2)])
        self.layer3 = DartsNode(prev_size=G2, init_tau=init_tau)

    def forward(self, X: torch.Tensor) -> torch.Tensor:
        L1 = torch.stack([n.forward(X) for n in self.layer1], dim=1)
        L2 = torch.stack([n.forward(L1) for n in self.layer2], dim=1)
        out = self.layer3.forward(L2)
        return out

    # --- discretisation ---
    def discretise(self, beta_threshold: float = 0.5) -> dict:
        """Return a dict describing the argmax architecture."""
        def describe(node: DartsNode):
            op = GATE_LIBRARY[int(torch.argmax(node.alpha).item())]
            sel = (torch.sigmoid(node.beta) > beta_threshold).nonzero(as_tuple=True)[0].tolist()
            gamma_pick = int(torch.argmax(node.gamma).item())
            return {"op": op, "selected": sel, "gamma_pick": gamma_pick,
                    "alpha": node.alpha.detach().cpu().numpy().tolist(),
                    "log_tau": float(node.log_tau.item())}

        return {
            "layer1": [describe(n) for n in self.layer1],
            "layer2": [describe(n) for n in self.layer2],
            "layer3": describe(self.layer3),
        }


def evaluate_discrete(circuit_desc: dict, X: np.ndarray, y: np.ndarray) -> dict:
    """Hard-evaluate the discretised circuit on data X (batch, N) -> y_hat in {0,1}."""
    def gate(op: str, sel_vals: np.ndarray, gamma_pick_val: float | None) -> float:
        # sel_vals: 1-D array of {0,1} selected input values
        if op == "AND":
            return 1.0 if (len(sel_vals) > 0 and np.all(sel_vals == 1)) else 0.0
        if op == "OR":
            return 1.0 if (len(sel_vals) > 0 and np.any(sel_vals == 1)) else 0.0
        if op == "XOR":
            return float(int(sel_vals.sum()) % 2) if len(sel_vals) > 0 else 0.0
        if op == "MAJ3":
            return 1.0 if int(sel_vals.sum()) >= 2 else 0.0
        if op == "MAJ5":
            return 1.0 if int(sel_vals.sum()) >= 3 else 0.0
        if op == "MAJ7":
            return 1.0 if int(sel_vals.sum()) >= 4 else 0.0
        if op == "ID":
            return float(gamma_pick_val) if gamma_pick_val is not None else 0.0
        if op == "NOT":
            return 1.0 - float(gamma_pick_val) if gamma_pick_val is not None else 1.0
        raise ValueError(f"unknown op {op}")

    n = X.shape[0]
    correct = 0
    preds = np.zeros(n, dtype=np.float32)
    for i in range(n):
        v0 = X[i].astype(np.float32)
        # Layer 1
        L1 = np.zeros(len(circuit_desc["layer1"]), dtype=np.float32)
        for j, nd in enumerate(circuit_desc["layer1"]):
            sel = nd["selected"]
            sel_vals = v0[sel] if len(sel) > 0 else np.array([], dtype=np.float32)
            gp = v0[nd["gamma_pick"]]
            L1[j] = gate(nd["op"], sel_vals, gp)
        # Layer 2
        L2 = np.zeros(len(circuit_desc["layer2"]), dtype=np.float32)
        for j, nd in enumerate(circuit_desc["layer2"]):
            sel = nd["selected"]
            sel_vals = L1[sel] if len(sel) > 0 else np.array([], dtype=np.float32)
            gp = L1[nd["gamma_pick"]]
            L2[j] = gate(nd["op"], sel_vals, gp)
        # Layer 3
        nd = circuit_desc["layer3"]
        sel = nd["selected"]
        sel_vals = L2[sel] if len(sel) > 0 else np.array([], dtype=np.float32)
        gp = L2[nd["gamma_pick"]]
        out = gate(nd["op"], sel_vals, gp)
        preds[i] = out
        if int(out) == int(y[i]):
            correct += 1
    return {"accuracy": correct / n, "preds": preds}


def train_one(target_fn, N: int, G1: int, G2: int, epochs: int, lr: float, seed: int,
              tau_anneal: bool = True, batch_size: int = 256, verbose: bool = False):
    torch.manual_seed(seed)
    np.random.seed(seed)

    n_all = np.arange(2 ** N, dtype=np.int64)
    y_all = target_fn(n_all)
    X_all = n_to_bits(n_all, N)
    X_t = torch.from_numpy(X_all)
    y_t = torch.from_numpy(y_all)

    circuit = DartsCircuit(N=N, G1=G1, G2=G2, init_tau=2.0)
    opt = torch.optim.Adam(circuit.parameters(), lr=lr)

    losses = []
    accs_train = []
    n_examples = X_t.shape[0]
    log_eps = 1e-6

    for ep in range(epochs):
        perm = torch.randperm(n_examples)
        ep_losses = []
        for start in range(0, n_examples, batch_size):
            idx = perm[start:start + batch_size]
            X_b = X_t[idx]
            y_b = y_t[idx]
            y_hat = circuit.forward(X_b)
            loss = F.binary_cross_entropy(y_hat.clamp(log_eps, 1 - log_eps), y_b)
            opt.zero_grad()
            loss.backward()
            opt.step()
            ep_losses.append(loss.item())
        epoch_loss = float(np.mean(ep_losses))
        losses.append(epoch_loss)

        # Anneal temperature: increase tau over training.
        if tau_anneal:
            with torch.no_grad():
                target_log_tau = float(np.log(2.0 + 8.0 * (ep + 1) / epochs))
                for module in circuit.modules():
                    if isinstance(module, DartsNode):
                        module.log_tau.data = module.log_tau.data * 0.99 + target_log_tau * 0.01

        if verbose and (ep + 1) % max(1, epochs // 10) == 0:
            with torch.no_grad():
                y_hat_full = circuit.forward(X_t).cpu().numpy()
                acc = float(np.mean((y_hat_full > 0.5).astype(np.float32) == y_all))
            accs_train.append(acc)
            print(f"  ep {ep+1:4d}/{epochs}  loss={epoch_loss:.4f}  acc(soft>0.5)={acc:.4f}")

    # Final continuous accuracy.
    circuit.eval()
    with torch.no_grad():
        y_hat_final = circuit.forward(X_t).cpu().numpy()
    soft_acc = float(np.mean((y_hat_final > 0.5).astype(np.float32) == y_all))

    desc = circuit.discretise(beta_threshold=0.5)
    discrete_eval = evaluate_discrete(desc, X_all, y_all)
    discrete_acc = discrete_eval["accuracy"]

    return {
        "losses": losses,
        "soft_acc": soft_acc,
        "discrete_acc": discrete_acc,
        "discrete_desc": desc,
        "discrete_preds_train": discrete_eval["preds"].tolist(),
        "y_train": y_all.tolist(),
        "soft_preds_train": y_hat_final.tolist(),
    }


def evaluate_extrapolation(desc: dict, n_lo: int, n_hi: int, N: int, target_fn) -> dict:
    n_arr = np.arange(n_lo, n_hi, dtype=np.int64)
    y_arr = target_fn(n_arr)
    # For inputs n with bit-length > N, we *truncate* to N low bits.
    X_arr = n_to_bits(n_arr & ((1 << N) - 1), N)
    res = evaluate_discrete(desc, X_arr, y_arr)
    return {
        "n_lo": int(n_lo),
        "n_hi": int(n_hi),
        "discrete_acc": float(res["accuracy"]),
        "y_density": float(np.mean(y_arr)),
    }


def random_bool_target(seed: int, N: int):
    """Return a callable n_array -> {0,1} that is a fixed random Boolean function on [0, 2^N)."""
    rng = np.random.RandomState(seed)
    table_lo = rng.rand(2 ** N) < 0.137  # PRIMES density at N=12 is 564/4096 ≈ 0.1377

    def fn(n_array: np.ndarray) -> np.ndarray:
        out = np.zeros_like(n_array, dtype=np.float32)
        mask = (n_array >= 0) & (n_array < 2 ** N)
        out[mask] = table_lo[n_array[mask]].astype(np.float32)
        # For n outside [0, 2^N): use 0 (consistent with prime extrapolation: the bits cycle)
        return out

    return fn


def run_experiment(args):
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    results = {"config": vars(args), "primes": [], "controls": [], "extrap": {}}

    # ---- PRIMES seeds ----
    for seed in range(args.n_seeds):
        t0 = time.time()
        print(f"PRIMES seed {seed}")
        res = train_one(
            target_fn=chi_P,
            N=args.N,
            G1=args.G1,
            G2=args.G2,
            epochs=args.epochs,
            lr=args.lr,
            seed=seed,
            verbose=args.verbose,
        )
        res["seed"] = seed
        res["wallclock_s"] = time.time() - t0
        results["primes"].append({
            "seed": seed,
            "final_loss": res["losses"][-1],
            "min_loss": min(res["losses"]),
            "soft_acc": res["soft_acc"],
            "discrete_acc": res["discrete_acc"],
            "discrete_desc": res["discrete_desc"],
            "wallclock_s": res["wallclock_s"],
            "losses": res["losses"],
        })
        print(f"  -> final_loss={res['losses'][-1]:.4f} soft_acc={res['soft_acc']:.4f} discrete_acc={res['discrete_acc']:.4f} ({res['wallclock_s']:.1f}s)")

    # ---- CONTROL seeds (matched-density random function) ----
    for seed in range(args.n_seeds):
        t0 = time.time()
        print(f"CONTROL seed {seed}")
        ctrl_fn = random_bool_target(seed=1000 + seed, N=args.N)
        res = train_one(
            target_fn=ctrl_fn,
            N=args.N,
            G1=args.G1,
            G2=args.G2,
            epochs=args.epochs,
            lr=args.lr,
            seed=seed + 100,
            verbose=args.verbose,
        )
        res["seed"] = seed
        res["wallclock_s"] = time.time() - t0
        results["controls"].append({
            "seed": seed,
            "final_loss": res["losses"][-1],
            "min_loss": min(res["losses"]),
            "soft_acc": res["soft_acc"],
            "discrete_acc": res["discrete_acc"],
            "discrete_desc": res["discrete_desc"],
            "wallclock_s": res["wallclock_s"],
            "losses": res["losses"],
        })
        print(f"  -> final_loss={res['losses'][-1]:.4f} soft_acc={res['soft_acc']:.4f} discrete_acc={res['discrete_acc']:.4f} ({res['wallclock_s']:.1f}s)")

    # ---- Extrapolation: best PRIMES seed (lowest final loss), evaluate on [2^N, 2^N+1000)
    best_p = min(results["primes"], key=lambda r: r["final_loss"])
    extr = evaluate_extrapolation(
        desc=best_p["discrete_desc"],
        n_lo=2 ** args.N, n_hi=2 ** args.N + 1000,
        N=args.N, target_fn=chi_P)
    results["extrap"]["primes_best"] = extr
    print(f"PRIMES best seed extrapolation [2^N, 2^N+1000) accuracy: {extr['discrete_acc']:.4f} (target density {extr['y_density']:.4f})")

    # save
    with open(out_dir / "results.json", "w") as f:
        json.dump(results, f, indent=2, default=lambda o: o.tolist() if isinstance(o, np.ndarray) else o)
    print(f"\nResults saved to {out_dir}/results.json")
    return results


def parse():
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=12)
    p.add_argument("--G1", type=int, default=12)
    p.add_argument("--G2", type=int, default=12)
    p.add_argument("--epochs", type=int, default=150)
    p.add_argument("--lr", type=float, default=1e-2)
    p.add_argument("--n_seeds", type=int, default=3)
    p.add_argument("--out-dir", type=str,
                   default="/apps/aplikacijos/prime-research/experiments/ml/darts_primes/run")
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()


if __name__ == "__main__":
    args = parse()
    run_experiment(args)
