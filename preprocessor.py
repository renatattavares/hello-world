# Run preprocessor

import numpy as np
from mspreprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as msh
import time
import pdb
import mspreprocessor.geoUtil.geoTools as gtool

from math import pi, sqrt
from pymoab import core, types, rng, topo_util

start = time.time()
M = msh("25.h5m", dim = 3)
end = time.time()
print("The preprocessing step lasted {0}s".format(end-start))
