import numpy as np

configfile: "config.yaml"

rule netzbooster:
    input:
        network = config["network"],
        variable= "data/variables.csv"
    output:
        network= "networks/new/postnetwork{flex_cost}_sn20.nc",
        P= "results/new/P_netzbooster{flex_cost}.csv",
        pos= "results/new/p_pos{flex_cost}.csv",
        neg= "results/new/p_neg{flex_cost}.csv"
    threads: 4
    resources: mem=10000
    script: "scripts/netzbooster_sim.py"

rule netzbooster_sclopf:
    input:
        network = config["network"],
        variable= "data/variables.csv"
    output:
        network= "networks/new/postnetwork_sclopf_{flex_cost}_sn20_{tatl_factor}.nc",
        P= "results/new/P_netzbooster_sclopf_{flex_cost}_{tatl_factor}.csv",
        pos= "results/new/p_pos_sclopf_{flex_cost}_{tatl_factor}.csv",
        neg= "results/new/p_neg_sclopf_{flex_cost}_{tatl_factor}.csv"
    threads: 4
    resources: mem=10000
    script: "scripts/netzbooster_sim.py"

def flex_pus():
    cf = config["flex_cost"]
    return np.round(np.arange(cf["start"], cf["stop"] + cf["step"], cf["step"]),2)

def tatl_pus():
    cf = config["tatl_factor"]
    return np.round(np.arange(cf["start"], cf["stop"] + cf["step"], cf["step"]),2)

rule all:
    input:
      booster_P       = expand("results/new/P_netzbooster{flex_cost}.csv", flex_cost=flex_pus()),
      booster_pos     = expand("results/new/p_pos{flex_cost}.csv", flex_cost=flex_pus()),
      booster_neg     = expand("results/new/p_neg{flex_cost}.csv", flex_cost=flex_pus()),
      booster_network = expand("networks/new/postnetwork{flex_cost}_sn20.nc", flex_cost=flex_pus())

rule all_sclopf:
    input:
      booster_P       = expand("results/new/P_netzbooster_sclopf_{flex_cost}_{tatl_factor}.csv", flex_cost=flex_pus(), tatl_factor=tatl_plus()),
      booster_pos     = expand("results/new/p_pos_sclopf_{flex_cost}_{tatl_factor}.csv", flex_cost=flex_pus(), tatl_factor=tatl_plus()),
      booster_neg     = expand("results/new/p_neg_sclopf_{flex_cost}_{tatl_factor}.csv", flex_cost=flex_pus(), tatl_factor=tatl_plus()),
      booster_network = expand("networks/new/postnetwork_sclopf_{flex_cost}_sn20_{tatl_factor}.nc", flex_cost=flex_pus(), tatl_factor=tatl_plus())
