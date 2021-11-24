import pypsa
from pyomo.environ import (ConcreteModel, Var, NonNegativeReals, Constraint,
                           Reals, Suffix, Binary, SolverFactory)
from pypsa.opt import (l_constraint, l_objective, LExpression, LConstraint)
import pandas as pd

import logging
logger = logging.getLogger(__name__)


import pypsa
from pyomo.environ import (ConcreteModel, Var, NonNegativeReals, Constraint,
                           Reals, Suffix, Binary, SolverFactory)
from pypsa.opt import (l_constraint, l_objective, LExpression, LConstraint)
import pandas as pd

from collections.abc import Iterable

from pypsa.descriptors import get_extendable_i, get_non_extendable_i
from pypsa.pf import calculate_PTDF, _as_snapshots
from pypsa.linopt import set_conref, write_constraint, get_var, linexpr


from contingency_tatl import network_sclopf_tatl




n = pypsa.Network(snakemake.input[0])

# n.calculate_dependent_values()


def get_branch_outages(n):
    """Creates list of considered outages"""
    outages = n.lines.index[n.lines.index.str.contains("_outage")].union(n.lines.index[~n.lines.index.str.contains("_outage") & ~(n.lines.index + "_outage").isin(n.lines.index)])
    return [("Line", o) for o in outages]





solver_options = {"name": 'gurobi', "threads": 4, "method": 2,
                   "crossover": 0, "BarConvTol": 1.e-4,
                   "FeasibilityTol": 1.e-5, "AggFill": 0,
                   "PreDual": 0, "GURO_PAR_BARDENSETHRESH": 200}

sclopf_kwargs = {
        "pyomo": False,
        "branch_outages": get_branch_outages(n),
        "solver_name": 'gurobi',
        "solver_options": solver_options,
        "formulation": "kirchhoff",
    }




snapshots = n.snapshots

tatl = float(snakemake.wildcards.tatl_factor)

network_sclopf_tatl(n,snapshots=snapshots, **sclopf_kwargs,tatl=tatl)



n.export_to_netcdf(snakemake.output[0])
# n.export_to_netcdf("postnetwork_full1_sn20_tatl.nc")





