from pkg_resources import get_distribution

try:
    __version__ = get_distribution("minospb").version
except:
    __version__ = "local"


__all__ = [
    "bayestyper",
    "graphtyper",
    "masking",
    "sample_dir",
    "tasks",
    "utils",
]

from minospb import *
