# Run preprocessor

import numpy as np
from meshHandle.multiscaleMesh import FineScaleMeshMS as msh
import time
import pdb
import geoUtil.geoTools as gtool

from math import pi, sqrt
from pymoab import core, types, rng, topo_util


M = msh("20.h5m", dim = 3)
