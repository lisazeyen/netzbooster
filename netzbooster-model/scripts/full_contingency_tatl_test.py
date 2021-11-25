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


solver_options = snakemake.config["solver"]


sclopf_kwargs = snakemake.config["sclopf"]
sclopf_kwargs["branch_outages"] = get_branch_outages(n)
sclopf_kwargs["solver_options"] = solver_options


snapshots = n.snapshots


tatl = float(snakemake.wildcards.tatlfactor)


network_sclopf_tatl(n, snapshots=snapshots, **sclopf_kwargs, tatl=tatl)


n.export_to_netcdf(snakemake.output[0])




