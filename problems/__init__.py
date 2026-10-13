"""Dynamic multi-objective test problems.

Two families are provided:

* ``FDA1``, ``FDA3``, ``FDA4``  -- Farina/Deb/Amato (2004), verified formula
  by formula against PlatEMO master.
* ``DMOP1``/``DMOP2``/``DMOP3`` -- project *variants* used in the v2.4
  experiments (see the module docstrings; DMOP3 is static by construction and
  serves as a static control).
* ``DMOP1A``/``DMOP2A``/``DMOP3A`` -- the **authoritative** dMOP1/dMOP2/dMOP3
  of Goh & Tan (2009), added in v2.5 for the supplementary experiment.
"""
from .base import DynamicProblem
from .fda1 import FDA1
from .fda3 import FDA3
from .fda4 import FDA4
from .dmop1 import DMOP1
from .dmop2 import DMOP2
from .dmop3 import DMOP3
from .dmop1a import DMOP1A
from .dmop2a import DMOP2A
from .dmop3a import DMOP3A
from .true_fronts import get_true_pf

__all__ = ["DynamicProblem", "FDA1", "FDA3", "FDA4", "DMOP1", "DMOP2",
           "DMOP3", "DMOP1A", "DMOP2A", "DMOP3A", "get_true_pf"]
