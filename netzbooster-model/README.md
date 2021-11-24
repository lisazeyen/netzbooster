# Workflow Overview

rules compute_approximate_C
  input: "prenetwork.nc"
  output: xyz
  script: noteboosk/contingency_c.py
  
  
# Find approximate buffer capacity factor
Here, we find buffer capacity factor based on our proposed approximate approach, c_a, for each random subset(cluster).

rules compute_approximate_c
  input: "networks/prenetwork.nc"
  output: c_a_cluster.csv
  script: scripts/find_c_a.py


# Find robust buffer capacity factor
Here, we find buffer capacity factor based on our proposed robust approach, c_r, for the network.

rules compute_robust_c
  input: "networks/prenetwork.nc"
  output: c_r.csv
  script: scripts/find_c_r.py


# Find robust buffer capacity factor
Here, we find buffer capacity factor based on our proposed line-specific approach, c_l, for the network.

rules compute_robust_l
  input: "networks/prenetwork.nc"
  output: c_l.csv
  script: scripts/find_c_l.py

# Heuristic Security-Constrained LOPF

This repository is useful to test the operational costs of
a network with heuristic contingency factors (cf. rule `solve_heuristic_contingency`) versus the
fully security-constrained case (cf. rule `solve_full_contingency`), as well as how
secure the heuristic systems for specific outages (cf. rule `check_outage_flow`).

## Workflow Structure

![workflow](workflow.png)

Documentation inside `scripts/*.py`.

## Configuration Options

The configuration is set in `config.yaml`.

- `network:` specifies the path to the already solved PyPSA-Eur network. Ideal would be a 256-node European network with investments optimized for 80% renewables, which has lots of flow, but also some conventional generators to set high marginal prices.
- `load_shedding:` adds load shedding generators to the network if not already existing to guarantee feasibility.
- `rolling_horizon:` runs full security-constraint LOPF in batches of snapshots rather than all at once.
- `group_size:` specifies the number of snapshots that form a batch.
- `s_max_pu:` specifies which heuristic contingency factors you want to sweep through.
- `solver:` includes the solver options and parameters.
 
## Executing on Cluster

Mount the LSDF, do `kubectl-login`, go to project directory and execute:

```sh
kubedo --inception -i pypsa-eur-gurobi -- snakemake -j 99 all
```

## Caveats

Make sure you do not have any lines with `num_parallel=0` and infinite reactance `x=np.inf`.

Make sure you do not have sub networks with a single bus.

## Analysis

You can get the operational costs of a network with

```py
import pypsa
n = pypsa.Network("results/postnetwork....nc")

(n.generators.marginal_cost * (n.snapshot_weightings @ n.generators_t.p)).sum()
```

You can check out the line loadings at the different outages with

```py
import pandas as pd

pd.read_csv("results/outage_line_loading_heur....csv", index_col=[0,1,2], parse_dates=True)
```

This file has multiple index levels: the first is the snapshot, the second the component type, the third the line index. The columns denote the outage of the line. The column `base` shows the flows under no outage conditions.

You can compare this data to the line capacities in `n.lines.s_nom` to find any overloadings.


## Plot



## TODOs

- [ ] Approximate, robust, and line-specific c-factors from the polygon analysis in https://git.scc.kit.edu/FN/contingency, but only once the group selection is improved.