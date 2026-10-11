"""Performance indicators and statistics for dynamic MOO experiments."""
from .igd import igd
from .hv import hv, get_reference_point, REF_POINTS
from .reference_front import reference_front
from .statistics import summarize, wilcoxon_holm
from .friedman_test import friedman_test, nemenyi_posthoc, friedman_nemenyi

__all__ = ["igd", "hv", "get_reference_point", "REF_POINTS",
           "reference_front", "summarize", "wilcoxon_holm",
           "friedman_test", "nemenyi_posthoc", "friedman_nemenyi"]
