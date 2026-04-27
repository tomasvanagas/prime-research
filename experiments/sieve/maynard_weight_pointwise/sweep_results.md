# sweep.py — results pointer

Driver script that runs `maynard_weight_pointwise.py` across the
(N, θ, F, H) configuration grid (92 configs total). Produces 92
`sweep_*.json` files plus a `sweep_summary.json` aggregate.

**Configuration grid:**
- `N ∈ {10^4, 10^4.5, 10^5}` × `θ ∈ {0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40}`
  × `F ∈ {selberg_gpy, linear, constant, maynard_sym}` (84 configs).
- Plus `H ∈ {(0,4,6), (0,2,6,8,12)}` × `θ ∈ {0.20, 0.30}` ×
  `N ∈ {10^4, 10^5}` (8 configs).

**Aggregated results and analysis:** see
`maynard_weight_pointwise_results.md` (this directory).

The sweep table itself is reproduced under §"Pointwise AUC" in the
main results file. Full per-config records are in `sweep_*.json` for
inspection / reanalysis.
