"""Preprocess network."""

__author__ = "Fabian Neumann (KIT)"
__copyright__ = f"Copyright 2020, {__author__}, GNU GPL 3"

import pypsa
import pandas as pd
import numpy as np

def add_NW_line():
    n.madd("Line",
           names=["NW"],
           num_parallel=2.67036,
           bus0=["DE0 11"],
           bus1=["DE0 15"],
           length=[147.677127],
           type=['Al/St 240/40 4-bundle 380.0'],
           carrier=['AC'],
           capital_cost = [15.35691],
           s_nom=[4534.545478]
    ) #same capacity/num_parallel as line "12", cost+length those of line "12" multiplied by 2

def prepare_network(n, Co2L_factor=1.):
    """
    Several preprocessing steps including load shedding and
    fixing optimised capacities.
    """

    n.lines.s_max_pu = 1.0

    # factor 10 because inital network is given with Co2L = 0.1.
    print(f"Setting a global Co2 emission target of {100*Co2L_factor} %")
    n.global_constraints.at["CO2Limit", "constant"] *= (10*Co2L_factor)
    
#     n.generators.loc[n.generators.carrier=='OCGT','p_nom_extendable'] = True
#     n.storage_units.p_nom_extendable = True
#     n.links.p_nom_extendable = True
#     n.lines.s_nom_extendable = True

#     n.generators.p_nom = n.generators.p_nom_opt
#     n.storage_units.p_nom = n.storage_units.p_nom_opt
#     n.links.p_nom = n.links.p_nom_opt
#     n.lines.s_nom = n.lines.s_nom_opt
    
#     n.remove("GlobalConstraint", "CO2Limit")
    
#     co2_price = 25
#     co2_pr = float(snakemake.wildcards.co2_price)
#     co2_costs=n.generators.carrier.map(n.carriers.co2_emissions)*co2_pr/n.generators.efficiency
    
#     n.generators.marginal_cost += co2_costs
    

    if snakemake.config["load_shedding"] and "load" not in n.carriers.index:
        n.add("Carrier", "load")
        n.madd(
            "Generator",
            n.buses.index,
            " load",
            bus=n.buses.index,
            carrier="load",
            sign=1e-3,  # measure p and p_nom in kW
            marginal_cost=1e2,  # Eur/kWh
            p_nom=1e9,  # kW
        )

    return n


def split_outage_lines(n):
    """
    Separate one parallel line in each corridor for which we consider the outage.
    Special handling for 220kV lines lifted to 380kV.
    """

    def create_outage_lines(n, condition, num_parallel):
        lines_out = n.lines.loc[condition].copy()
        lines_out.s_nom = lines_out.s_nom * num_parallel / lines_out.num_parallel
        lines_out.num_parallel = num_parallel
        lines_out.index = [f"{i}_outage" for i in lines_out.index]
        return lines_out

    def adjust_nonoutage_lines(n, condition, num_parallel):
        nump = n.lines.loc[condition, "num_parallel"]
        n.lines.loc[condition, "s_nom"] *= (nump - num_parallel) / nump
        n.lines.loc[condition, "num_parallel"] -= num_parallel

    cond_220 = (n.lines.num_parallel > 0.5) & (n.lines.num_parallel < 1)
    cond_380 = n.lines.num_parallel > 1

    lines_out_220 = create_outage_lines(n, cond_220, 1 / 3)
    lines_out_380 = create_outage_lines(n, cond_380, 1.0)

    adjust_nonoutage_lines(n, cond_220, 1 / 3)
    adjust_nonoutage_lines(n, cond_380, 1)

    n.lines = pd.concat([n.lines, lines_out_220, lines_out_380])

    n.calculate_dependent_values()

    return n



def get_representative_snapshots(n, n_snapshots):
    assert n.storage_units.empty & n.stores.empty, "Stores and storage_units must be empty!"
    try:
        import tsam.timeseriesaggregation as tsam
    except:
        raise ModuleNotFoundError("Optional dependency 'tsam' not found."
                                  "Install via 'pip install tsam'")

    p_max_pu_norm = n.generators_t.p_max_pu.max()
    p_max_pu = n.generators_t.p_max_pu / p_max_pu_norm
    
    load_norm = n.loads_t.p_set.max()
    load = n.loads_t.p_set / load_norm

    raw = pd.concat([p_max_pu, load], axis=1, sort=False)

    solver_name = 'gurobi'

    agg = tsam.TimeSeriesAggregation(raw,
                                     hoursPerPeriod=1,
                                     noTypicalPeriods=n_snapshots,
                                     extremePeriodMethod='replace_cluster_center',
                                     solver=solver_name)

    clustered = agg.createTypicalPeriods()

    
    clustered_ind = raw.iloc[agg.clusterCenterIndices].index
    n.set_snapshots(clustered_ind)
    # this can be done nicer, probably, but does the job (for now):
    n.snapshot_weightings.objective = list(agg._clusterPeriodNoOccur.values())
    n.snapshot_weightings.stores = list(agg._clusterPeriodNoOccur.values())
    n.snapshot_weightings.generators = list(agg._clusterPeriodNoOccur.values())

    return n


# def apply_hacks(n):
#     """
#     Here's some space to do some hacking to a specific network.
#     """

#     to_remove = n.lines.loc[n.lines.s_nom_opt < 10].index
#     print("Removing following lines with small capacity:",to_remove)
#     n.mremove("Line", to_remove)

#     to_remove = n.lines.index[n.lines.x == np.inf]
#     print("Removing following lines with infinite reactance:",to_remove)
#     n.mremove("Line", to_remove)

#     bus_to_remove = "DE0 62"
#     n.remove("Bus", bus_to_remove)
#     n.mremove("StorageUnit", n.storage_units.loc[n.storage_units.bus == bus_to_remove].index)
#     n.mremove("Generator", n.generators.loc[n.generators.bus == bus_to_remove].index)
#     n.mremove("Load", n.loads.loc[n.loads.bus == bus_to_remove].index)

#     return n


if __name__ == "__main__":

    n = pypsa.Network(snakemake.input[0])

    add_NW_line()

    n = prepare_network(n, Co2L_factor = float(snakemake.wildcards["Co2L"].split('L')[1]))

#    n = apply_hacks(n)

    n = split_outage_lines(n)

    to_remove = n.lines.loc[n.lines.num_parallel==0].index
    print("Removing following lines with num_parallel = 0:",to_remove)
    n.mremove("Line", to_remove)
    
    n.mremove("StorageUnit", n.storage_units.index) 
    
    n.determine_network_topology()

    N_SNAPSHOTS = int(snakemake.wildcards["snapshots"].split('n')[1])
    n = get_representative_snapshots(n, N_SNAPSHOTS)


    n.export_to_netcdf(snakemake.output[0])
