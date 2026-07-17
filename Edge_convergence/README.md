# Edge_convergence

This directory contains edge-level convergence data for the initial retrospective free energy perturbation (FEP) calculations described in the manuscript.

The last-quarter median slope (LQMS) metric was calculated as used in Sampson et al. *JMB* (2024): it is the Theil-Sen slope over the last quarter of the trajectory, with a minimum window of 5 ns.
Because the FEP calculations were 10 ns long, the LQMS values reported here were calculated over the final 5 ns of simulation.

## Files

- `supplemental_edge_convergence_data.csv`: convergence metrics and threshold flags for each FEP edge.

## Columns

- `run_group`: calculation set represented by the row.
- `fmp_id`: identifier for the source FEP map.
- `edge_id`: hash-based identifier for the FEP edge.
- `n0`, `n1`: the two endpoint states defining the edge.
- `mutation`: mutation represented by the edge.
- `cpx_lqms_final`, `sol_lqms_final`, `ddg_lqms_final`: final LQMS values for the complex leg, solvent leg, and Bennett relative binding free energy, respectively.
- `cpx_lqvar_final`, `sol_lqvar_final`, `ddg_lqvar_final`: final last-quarter variance values for the complex leg, solvent leg, and Bennett relative binding free energy, respectively.
- `cpx_above_<threshold>`, `sol_above_<threshold>`, `ddg_above_<threshold>`: Boolean flags indicating whether the corresponding LQMS magnitude exceeds the named threshold (0.05, 0.1, 0.2, or 0.5).
- `any_above_<threshold>`: Boolean flags indicating whether the LQMS magnitude for any of the complex, solvent, or relative binding free-energy values exceeds the named threshold.

Missing values are represented by `NA`, and were due to missing automated analysis data from the FEP edge.
