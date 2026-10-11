"""Dynamic multi-objective test problems (FDA1, FDA3, FDA4, dMOP1-3)."""
from .base import DynamicProblem
from .fda1 import FDA1
from .fda3 import FDA3
from .fda4 import FDA4
from .dmop1 import DMOP1
from .dmop2 import DMOP2
from .dmop3 import DMOP3
from .true_fronts import get_true_pf

__all__ = ["DynamicProblem", "FDA1", "FDA3", "FDA4", "DMOP1", "DMOP2",
           "DMOP3", "get_true_pf"]
