import pypsa
from pypsa.linopt import (get_var, linexpr, define_constraints, write_objective,
                          define_variables, get_con, run_and_read_gurobi, join_exprs)
from pypsa.linopf import prepare_lopf
from pypsa.descriptors import expand_series, additional_linkports
from pypsa.pf import get_switchable_as_dense as get_as_dense, _as_snapshots
import pandas as pd
from numpy import inf
import numpy as np
import os
import logging
logger = logging.getLogger(__name__)
from tempfile import mkstemp
import gc
from vresutils import Dict
import yaml

import sys
sys.path.append('../netzbooster-model/scripts/')
from contingency_tatl import add_contingency_constraints_lowmem

if "snakemake" not in globals():
    from vresutils import Dict
    import yaml
    snakemake = Dict()
    snakemake = Dict()
    snakemake.config = Dict()
    with open('../../netzbooster/netzbooster-model/config.yaml', encoding='utf8') as f:
        snakemake.config = yaml.safe_load(f)
        config  = snakemake.config
        snakemake.input = Dict(
        network  = "../netzbooster-model/networks/prenetwork_Co2L0.0_sn3.nc",
        variable = "data/variables.csv"
        )
        snakemake.output = Dict(
            network = "results/postnetwork_sclopf_Co2L0.0_sn3"+
                      "_cost{flex_cost}_tatl1.0.nc",
            P       = "results/P_netzbooster_sclopf_Co2L0.0_sn3"+
                      "_cost10_tatl1.0.csv",
            pos     = "results/p_pos_sclopf_Co2L0.0_sn3"+
                       "_cost10_tatl1.0.csv",
            neg     = "results/p_neg_sclopf_Co2L0.0_sn3"+
                      "_cost10_tatl1.0.csv"
                    )
        snakemake.wildcards = Dict(
            tatlfactor = "1.0",
            flex_cost="10"
            )

#%%
# functions ------------------------------------------------------------------
def get_branch_outages(n):
    """Creates list of considered outages"""
    outages = n.lines.index[n.lines.index.str.contains("_outage")].union(n.lines.index[~n.lines.index.str.contains("_outage") & ~(n.lines.index + "_outage").isin(n.lines.index)])
    return [("Line", o) for o in outages]

def assign_solution_netzbooster(n, sns, variables_sol, constraints_dual,
                    keep_references=False, keep_shadowprices=None):
    """
    Helper function. Assigns the solution of a succesful optimization to the
    network.

    """

    def set_from_frame(pnl, attr, df):
        if attr not in pnl: #use this for subnetworks_t
            pnl[attr] = df.reindex(n.snapshots)
        elif pnl[attr].empty:
            pnl[attr] = df.reindex(n.snapshots)
        else:
            pnl[attr].loc[sns, :] = df.reindex(columns=pnl[attr].columns)

    pop = not keep_references

    def map_solution(c, attr):
        variables = get_var(n, c, attr, pop=pop)
        predefined = True
        if (c, attr) not in lookup.index:
            predefined = False
            n.sols[c] = n.sols[c] if c in n.sols else Dict(df=pd.DataFrame(), pnl={})
        n.solutions.at[(c, attr), 'in_comp'] = predefined
        if isinstance(variables, pd.DataFrame):
            # case that variables are timedependent
            n.solutions.at[(c, attr), 'pnl'] = True
            pnl = n.pnl(c) if predefined else n.sols[c].pnl
            values = variables.apply(lambda x: x.map(variables_sol))
            # values = variables.stack().map(variables_sol).unstack()
            if c in n.passive_branch_components and attr == "s":
                set_from_frame(pnl, 'p0', values)
                set_from_frame(pnl, 'p1', - values)
            elif c == 'Link' and attr == "p":
                set_from_frame(pnl, 'p0', values)
                for i in ['1'] + additional_linkports(n):
                    i_eff = '' if i == '1' else i
                    eff = get_as_dense(n, 'Link', f'efficiency{i_eff}', sns)
                    set_from_frame(pnl, f'p{i}', - values * eff)
                    pnl[f'p{i}'].loc[sns, n.links.index[n.links[f'bus{i}'] == ""]] = \
                        n.component_attrs['Link'].loc[f'p{i}','default']
            else:
                set_from_frame(pnl, attr, values)
        else:
            # case that variables are static
            n.solutions.at[(c, attr), 'pnl'] = False
            sol = variables.map(variables_sol)
            if predefined:
                non_ext = n.df(c)[attr]
                n.df(c)[attr + '_opt'] = sol.reindex(non_ext.index).fillna(non_ext)
            else:
                n.sols[c].df[attr] = sol

    n.sols = Dict()
    n.solutions = pd.DataFrame(index=n.variables.index, columns=['in_comp', 'pnl'])
    for c, attr in n.variables.index:
        map_solution(c, attr)

    # if nominal capcity was no variable set optimal value to nominal
    for c, attr in lookup.query('nominal').index.difference(n.variables.index):
        n.df(c)[attr+'_opt'] = n.df(c)[attr]

    # recalculate storageunit net dispatch
    if not n.df('StorageUnit').empty:
        c = 'StorageUnit'
        n.pnl(c)['p'] = n.pnl(c)['p_dispatch'] - n.pnl(c)['p_store']

    # duals
    if keep_shadowprices == False:
        keep_shadowprices = []

    sp = n.constraints.index
    if isinstance(keep_shadowprices, list):
        sp = sp[sp.isin(keep_shadowprices, level=0)]

    def map_dual(c, attr):
        # If c is a pypsa component name the dual is store at n.pnl(c)
        # or n.df(c). For the second case the index of the constraints have to
        # be a subset of n.df(c).index otherwise the dual is stored at
        # n.duals[c].df
        constraints = get_con(n, c, attr, pop=pop)
        is_pnl = isinstance(constraints, pd.DataFrame)
        sign = 1 if 'upper' in attr or attr == 'marginal_price' else -1
        n.dualvalues.at[(c, attr), 'pnl'] = is_pnl
        to_component = c in n.all_components
        if is_pnl:
            n.dualvalues.at[(c, attr), 'in_comp'] = to_component
            # changed for netzbooster
            duals = constraints.apply(lambda x: x.map(sign * constraints_dual))
            if c not in n.duals and not to_component:
                n.duals[c] = Dict(df=pd.DataFrame(), pnl={})
            pnl = n.pnl(c) if to_component else n.duals[c].pnl
            set_from_frame(pnl, attr, duals)
        else:
            # here to_component can change
            duals = constraints.map(sign * constraints_dual)
            if to_component:
                to_component = (duals.index.isin(n.df(c).index).all())
            n.dualvalues.at[(c, attr), 'in_comp'] = to_component
            if c not in n.duals and not to_component:
                n.duals[c] = Dict(df=pd.DataFrame(), pnl={})
            df = n.df(c) if to_component else n.duals[c].df
            df[attr] = duals

    n.duals = Dict()
    n.dualvalues = pd.DataFrame(index=sp, columns=['in_comp', 'pnl'])
    # extract shadow prices attached to components
    for c, attr in sp:
        map_dual(c, attr)

    #correct prices for snapshot weightings
    n.buses_t.marginal_price.loc[sns] = n.buses_t.marginal_price.loc[sns].divide(n.snapshot_weightings.objective.loc[sns],axis=0)

    # discard remaining if wanted
    if not keep_references:
        for c, attr in n.constraints.index.difference(sp):
            get_con(n, c, attr, pop)

    #load
    if len(n.loads):
        set_from_frame(n.pnl('Load'), 'p', get_as_dense(n, 'Load', 'p_set', sns))

    #clean up vars and cons
    for c in list(n.vars):
        if n.vars[c].df.empty and n.vars[c].pnl == {}: n.vars.pop(c)
    for c in list(n.cons):
        if n.cons[c].df.empty and n.cons[c].pnl == {}: n.cons.pop(c)

    # recalculate injection
    ca = [('Generator', 'p', 'bus' ), ('Store', 'p', 'bus'),
          ('Load', 'p', 'bus'), ('StorageUnit', 'p', 'bus'),
          ('Link', 'p0', 'bus0'), ('Link', 'p1', 'bus1')]
    for i in additional_linkports(n):
        ca.append(('Link', f'p{i}', f'bus{i}'))

    sign = lambda c: n.df(c).sign if 'sign' in n.df(c) else -1 #sign for 'Link'
    n.buses_t.p = pd.concat(
            [n.pnl(c)[attr].mul(sign(c)).rename(columns=n.df(c)[group])
             for c, attr, group in ca], axis=1).groupby(level=0, axis=1).sum()\
            .reindex(columns=n.buses.index, fill_value=0)

    def v_ang_for_(sub):
        buses_i = sub.buses_o
        if len(buses_i) == 1:
            return pd.DataFrame(0, index=sns, columns=buses_i)
        sub.calculate_B_H(skip_pre=True)
        Z = pd.DataFrame(np.linalg.pinv((sub.B).todense()), buses_i, buses_i)
        Z -= Z[sub.slack_bus]
        return n.buses_t.p.reindex(columns=buses_i) @ Z
    n.buses_t.v_ang = (pd.concat([v_ang_for_(sub) for sub in n.sub_networks.obj],
                                  axis=1)
                      .reindex(columns=n.buses.index, fill_value=0))

############################################################################
# ------------ Netzbooster modifications -----------------------------------
def add_to_objective(n, snapshots):
    """
    adds to the pypsa objective function the costs for the Netzbooster,
    these consists of:
        (a) costs for the capacities of the Netzbooster
        (b) costs for compensation dispatch (e.g. paid to DSM consumers, for storage)

    i => bus
    t => snapshot
    k => outage
    """

    logger.info("define objective")

    n.determine_network_topology()
    sub = n.sub_networks.at["0", "obj"]
    sub.calculate_PTDF()
    sub.calculate_BODF()

    buses = sub.buses_i()
    lines = sub.lines_i()
    outages = [l for l in lines if "outage" in l]

    coefficient1 = float(snakemake.wildcards.flex_cost) #18000
    coefficient2 = 10 #float(snakemake.wildcards.flex_cost)*0.001 #1

    #  (a) costs for Netzbooster capacities
    # sum_i ( c_pos * P(i)_pos + c_neg * P(i)_neg)
    netzbooster_cap_cost = linexpr(
                                (coefficient1, get_var(n, "Bus", "P_pos").loc[buses]),
                                (coefficient1, get_var(n, "Bus", "P_neg").loc[buses])
                                ).sum()

    write_objective(n, netzbooster_cap_cost)


    # (b) costs for compensation dispatch (paid to DSM consumers/running costs for storage)
    # sum_(i, t, k) ( o_pos * p(i, t, k)_pos + o_neg * p(i, t, k)_pos)

    compensate_p_cost = linexpr(
                                (coefficient2, get_var(n, "Bus", "p_pos").loc[snapshots, buses]),
                                (coefficient2, get_var(n, "Bus", "p_neg").loc[snapshots, buses]),
                                ).sum().sum()

    write_objective(n, compensate_p_cost)



def netzbooster_constraints(n, snapshots):
    """
    add to the LOPF the additional Netzbooster constraints
    """

    # Define indices & parameters
    n.determine_network_topology()
    sub = n.sub_networks.at["0", "obj"]
    sub.calculate_PTDF()
    sub.calculate_BODF()

    buses = sub.buses_i()
    lines = sub.lines_i()
    outages = [l for l in lines if "outage" in l]

    ptdf = pd.DataFrame(sub.PTDF, index=lines, columns=buses)
    lodf = pd.DataFrame(sub.BODF, index=lines, columns=lines)

    logger.info("define variables")

    # i = buses
    # t = snapshots
    # k = outages

    # P_i + >= 0
    define_variables(n, 0, inf, "Bus", "P_pos", axes=[buses])
    # P_i - >= 0
    define_variables(n, 0, inf, "Bus", "P_neg", axes=[buses])
    # p_{i,t,k} + >= 0
    define_variables(n, 0, inf, "Bus", "p_pos", axes=[snapshots,
                                                      pd.MultiIndex.from_product([buses, outages])])
    # p_{i,t,k} + >= 0
    define_variables(n, 0, inf, "Bus", "p_neg", axes=[snapshots,
                                                      pd.MultiIndex.from_product([buses, outages])])

    logger.info("define constraints")

    # # Define constraints

    # ########### p_pos(i,k,t) <= P_pos(i)  ##########
    P_pos_extend = expand_series(get_var(n, "Bus", "P_pos").loc[buses],
                                 snapshots).T.reindex(
                                     pd.MultiIndex.from_product([buses, outages]), axis=1, level=0)
    lhs_1 = linexpr(
                    (1,  P_pos_extend[buses]),
                    (-1, get_var(n, "Bus", "p_pos").loc[snapshots, buses])
                   )
    define_constraints(n, lhs_1, ">=", 0., "Bus", "UpLimitPos")


    # ######### p_neg(i,k,t) <= P_neg(i)  #########
    P_neg_extend = expand_series(get_var(n, "Bus", "P_neg").loc[buses],
                                 snapshots).T.reindex(
                                     pd.MultiIndex.from_product([buses, outages]), axis=1, level=0)
    lhs_2 = linexpr(
                    (1,  P_neg_extend[buses]),
                    (-1, get_var(n, "Bus", "p_neg").loc[snapshots, buses])
                   )


    define_constraints(n, lhs_2, ">=", 0., "Bus", "UpLimitNeg")

    # ######## sum(i, p_pos(i,k,t) - p_neg(i,k,t)) = 0  #########
    lhs_3 = linexpr(
                    (1,  get_var(n, "Bus", "p_pos")),
                    (-1,  get_var(n, "Bus", "p_neg"))
                   ).groupby(level=1, axis=1).sum()
    define_constraints(n, lhs_3, "<=", 0.1, "Bus", "FlexBal_ub")
    define_constraints(n, lhs_3, ">=", -0.1, "Bus", "FlexBal_lb")

    # ## |f(l,t) + LODF(l,k) * f(k,t) + sum_i ((PTDF(l,i) +  LODF(l,k) * PTDF(k,i)) * (p_neg(i,k, t) - p_pos(i,k, t)))| <= F(l)  ########
    # f(l,t) + LODF(l,k) * f(k,t) + sum_i ((PTDF(l,i) +  LODF(l,k) * PTDF(k,i))* (p_neg(i,k, t) - p_pos(i,k, t))) <= F(l)

    # f(l,t) here: f pandas DataFrame(index=snapshots, columns=lines)
    f = get_var(n,"Line", "s")
    # F(l) -> line capacity for n.lines.s_nom_extendable=False
    F = n.lines.s_nom
    # p pos (index=snapshots, columns=[bus, outage])
    p_pos = get_var(n, "Bus", "p_pos")
    # p neg (index=snapshots, columns=[bus, outage])
    p_neg = get_var(n, "Bus", "p_neg")


    for k in outages:
        # a = (PTDF(l,i) +  LODF(l,k) * PTDF(k,i))
        a = ptdf + expand_series(lodf[k], buses).mul(ptdf.loc[k])
        # expand a with snapshots [index=snapshots, columns=[bus, lines]]
        a_t = expand_series(a.T.stack(), snapshots).T
        # all lines except outage linke k
        l = lines.difference(pd.Index([k]))
        # get PTDF without outage line k (index=snapshots, columns=[buses, lines])
        a_no_k = a_t.reindex(l, axis=1, level=1)
        # expand lodf for linear expression
        lodf_t = expand_series(lodf.loc[l,k], snapshots).T
        # expand p_neg with lines (index=snapshots, columns=[buses, lines])
        p_neg_l = p_neg.xs(k, level=1, axis=1).reindex(a_no_k.columns, axis=1, level=0)
        # expand p_pos with lines (index=snapshots, columns=[buses, lines])
        p_pos_l = p_pos.xs(k, level=1, axis=1).reindex(a_no_k.columns, axis=1, level=0)
        # expand flow on outage line
        f_k = expand_series(f[k], l)

        # sum_i (a * (p_pos(i,k, t) - p_neg(i,k, t)))
        lhs = linexpr((a_no_k, p_pos_l),    # sum_i (a(l,i) * (p_pos(i,k, t))
                    (-1 * a_no_k,  p_neg_l)    # sum_i (-a(l,i) * (p_neg(i,k, t))
                    ).groupby(level=1, axis=1).sum()

        # f(l,t) + LODF(l,k) * f(k,t)
        lhs += linexpr((1, f[l]),             # f(l,t)
                       (lodf_t, f_k))   # LODF(l,k) * f(k,t)


        rhs = expand_series(F.loc[l], snapshots).T

        define_constraints(n, lhs, "<=", rhs, "booster1", "capacity_limit_lb_{}".format(k))

        define_constraints(n, lhs, ">=", -rhs, "booster2", "capacity_limit_ub_{}".format(k))



def extra_functionality(n,snapshots):

    netzbooster_constraints(n, snapshots)

    logger.info("Add contingency constraints")
    add_contingency_constraints_lowmem(n, snapshots,tatl)

    add_to_objective(n, snapshots)



# ---------------- LOPF with netzbooster -------------------------------------
def netzbooster_sclopf(n, snapshots, extra_functionality,
                     pyomo=False, formulation="kirchhoff", keep_references=True,
                     skip_objective=False, keep_files=False, solver_dir=None,
                     solver_logfile=None,
                     keep_shadowprices=['Bus', 'Line', 'Transformer', 'Link', 'GlobalConstraint'],
                     store_basis=False, warmstart=False, branch_outages=None, tatl=None):

    # prepare network
    snapshots = _as_snapshots(n, snapshots)
    n.calculate_dependent_values()
    n.determine_network_topology()

    # for conitengency constraints
    passive_branches = n.passive_branches()

    if branch_outages is None:
        branch_outages = passive_branches.index

    # save to network for extra_functionality
    n._branch_outages = branch_outages

    # to avoid numerical problems with small values
    for df in (n.generators_t.p_max_pu, n.storage_units_t.inflow):
            df.where(df>snakemake.config["solver_options"]['clip_p_max_pu'],
                     other=0., inplace=True)

    # formulate optimisation problem ---------------------------------------------
    # include standard pypsa constraints and objective + netzbooster (extra_functionality)
    logger.info("Prepare linear problem")
    # careful problem and constraints are stored in /tmp in a .txt file
    # with only 2 snapshots this file has already ~100 MB
    fdp, problem_fn = prepare_lopf(n, snapshots, keep_files, skip_objective,
                                   extra_functionality, solver_dir)
    fds, solution_fn = mkstemp(prefix='pypsa-solve', suffix='.sol', dir=solver_dir)

    # solve the linear problem --------------------------------------------------
    logger.info(f"Solve linear problem using {solver_name.title()} solver")
    solve = eval(f'run_and_read_{solver_name}')
    # with 2 snapshots this takes about 4 minutes (230.99 seconds)
    res = solve(n, problem_fn, solution_fn, solver_logfile,
                solver_options, warmstart, store_basis)
    status, termination_condition, variables_sol, constraints_dual, obj = res

    # close the LP files for solving
    if not keep_files:
        os.close(fdp); os.remove(problem_fn)
        os.close(fds); os.remove(solution_fn)

    # save objective value in network
    n.objective = obj
    assign_solution_netzbooster(n, snapshots, variables_sol, constraints_dual,
                    keep_references=keep_references,
                    keep_shadowprices=keep_shadowprices)
    gc.collect()

    return n

# ------ save results --------------------------------------------------------
def save_results_as_csv(n, snapshots):

    # Netzbooster capacity P_pos (pandas Series, index = Buses)
    P_pos = n.sols.Bus.df["P_pos"]
    # Netzbooster capacity P_neg (pandas Series, index = Buses)
    P_neg = n.sols.Bus.df["P_neg"]
    # pd.concat([P_pos, P_neg],axis=1).to_csv("results/P_netzbooster_capacity.csv")
    pd.concat([P_pos, P_neg],axis=1).to_csv(snakemake.output.P)

    # Netzbooster shadow price p_pos (only defined for considered snapshots, others are NaN values)
    # Dataframe(index=snapshots, columns=[Bus, outage])
    p_pos = n.sols.Bus.pnl["p_pos"].loc[snapshots]
    # p_pos.to_csv("results/p_pos.csv")
    p_pos.to_csv(snakemake.output.pos)
    p_neg = n.sols.Bus.pnl["p_neg"].loc[snapshots]
    # p_neg.to_csv("results/p_neg.csv")
    p_neg.to_csv(snakemake.output.neg)

    # save constraints with multiindex as csv and remove from network
    keys_to_remove = []
    for key in n.buses_t.keys():
        # n.buses_t[key].loc[snapshots].to_csv("results/{}.csv".format(key))
        if isinstance(n.buses_t[key].columns, pd.MultiIndex):
            keys_to_remove.append(key)

    # remove dataframes with multi-columns to avoid errors when exporting the network
    for key in keys_to_remove:
        n.buses_t.pop(key, None)


# import network
n = pypsa.Network(snakemake.input.network)

n.lines.s_max_pu = 1

lookup = pd.read_csv(snakemake.input.variable,
                     index_col=['component', 'variable'])

solver_options = snakemake.config["solver"]

# solver_name = "gurobi"
solver_name = solver_options.pop("name")

# for testing only consider 2 snapshots
snapshots = n.snapshots #[:15]

# branch outages
branch_outages=get_branch_outages(n)
tatl = float(snakemake.wildcards.tatlfactor)

#fix n._multi_invest problem by setting it manually to 0
n._multi_invest = 0

# run lopf with netzbooster constraints and modified objective
n = netzbooster_sclopf(n, snapshots, extra_functionality,
                       branch_outages=branch_outages, tatl=tatl)

# Netzbooster capacity P_pos (pandas Series, index = Buses)
P_pos = n.sols.Bus.df["P_pos"]
# Netzbooster capacity P_neg (pandas Series, index = Buses)
P_neg = n.sols.Bus.df["P_neg"]
# small p are saved here (pandas Dataframe, index=snapshots, columns=MultiIndex([bus, outage]))
p_pos = n.sols.Bus.pnl["p_pos"].loc[snapshots]
p_neg = n.sols.Bus.pnl["p_neg"].loc[snapshots]

save_results_as_csv(n, snapshots)

# n.export_to_netcdf("networks/postnetwork4_sn50.nc")
n.export_to_netcdf(snakemake.output.network)
#%% check constraints

#(1) P_pos >= 0
if not all(P_pos>=0):
    print("Warning not all capital P_pos>=0")
#(2) P_neg >= 0
if not all(P_neg>=0):
    print("Warning not all capital P_neg>=0")
#(3) 0<=p_pos<=P_pos
if not all((p_pos>=0) & (p_pos <= expand_series(P_pos.reindex(p_pos.columns, level=0),snapshots).T)):
    print("Warning not all small p_pos in 0<=p_pos<=P_pos")
#(4) 0<=p_neg<= P_neg
if not all((p_neg>=0) & (p_neg <= expand_series(P_neg.reindex(p_neg.columns, level=0),snapshots).T)):
    print("Warning not all small p_neg in 0<=p_neg<=P_neg")

#(5) sum_i (p_pos(i,t,k) - p_neg(i,t,k)) = 0
if not all(round((p_pos-p_neg).groupby(level=1, axis=1).sum(), ndigits=6)==0):
    print("constraint sum_i (p_pos(i,t,k) - p_neg(i,t,k)) = 0 is violated")

#(6) |f(l,t) + LODF(l,k) * f(k,t) + sum_i ((PTDF(l,i) +  LODF(l,k) * PTDF(k,i)) * (p_neg(i,k, t) - p_pos(i,k, t)))| <= F(l)  ########
F = n.lines.s_nom
f = n.lines_t.p0
n.determine_network_topology()
sub = n.sub_networks.at["0", "obj"]
sub.calculate_PTDF()
sub.calculate_BODF()

buses = sub.buses_i()
lines = sub.lines_i()
outages = [l for l in lines if "outage" in l]

ptdf = pd.DataFrame(sub.PTDF, index=lines, columns=buses)
lodf = pd.DataFrame(sub.BODF, index=lines, columns=lines)

ptdf_t = expand_series(ptdf.T.stack(), snapshots).T

for k in outages:

    # get values for outage k ----------------------------------------------
    # all lines except outage linke k
    l = lines.difference(pd.Index([k]))
    # a = (PTDF(l,i) +  LODF(l,k) * PTDF(k,i))
    a = ptdf + expand_series(lodf[k], buses).mul(ptdf.loc[k])
    # expand a with snapshots [index=snapshots, columns=[bus, lines]]
    a_t = expand_series(a.T.stack(), snapshots).T
    # all lines except outage linke k
    l = lines.difference(pd.Index([k]))
    # get PTDF without outage line k (index=snapshots, columns=[buses, lines])
    a_no_k = a_t.reindex(l, axis=1, level=1)
    # expand lodf for linear expression
    lodf_t = expand_series(lodf.loc[l,k], snapshots).T
    # expand p_neg with lines (index=snapshots, columns=[buses, lines])
    p_neg_l = p_neg.xs(k, level=1, axis=1).reindex(a_no_k.columns, axis=1, level=0)
    # expand p_pos with lines (index=snapshots, columns=[buses, lines])
    p_pos_l = p_pos.xs(k, level=1, axis=1).reindex(a_no_k.columns, axis=1, level=0)
    # expand flow on outage line
    f_k = expand_series(f[k], l)

    # terms of the final constraint ----------------------------------------
    # first term f(l,t) + LODF(l,k) * f(k,t)
    term1 = f.loc[snapshots, l] + lodf_t*f_k.loc[snapshots]
    # second term sumation over buses
    # sum_i (a (l,i) * (p_neg(i,k, t) - p_pos(i,k, t)))
    term2 = (a_no_k * (p_pos_l-p_neg_l)).groupby(level=1, axis=1).sum()

    if not all(abs(term1 + term2) <= expand_series(F.loc[l], snapshots).T):
        print("Warning constraint |f(l,t) + LODF(l,k) * f(k,t) + sum_i (PTDF(l,i) * (p_neg(i,k, t) - p_pos(i,k, t)))| <= F(l) "
              "is violated for outage {}".format(k))
