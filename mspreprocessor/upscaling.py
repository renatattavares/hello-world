# UPSCALLING OF STRUCTURED MESHES
import numpy as np
import time
import pdb
import geoUtil.geoTools as gtool
import xlsxwriter
from math import pi, sqrt
from pymoab import rng, types
from meshHandle.multiscaleMesh import FineScaleMeshMS as msh
from meshHandle.corePymoab import CoreMoab as core
from tpfa.boundary_conditions import BoundaryConditions
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.linalg import spsolve

def equiv_perm(k1, k2):
    return (2*k1*k2)/(k1 + k2)

def centroid_dist(c1, c2):
    return ((c1-c2)**2).sum()

print("Initializating mesh")
start = time.time()
dx, dy, dz = 1, 1, 1
nx, ny, nz = 10, 10, 10
num_elements = nx*ny*nz
M = msh("malha-teste.h5m", dim = 3)
vec = np.arange(len(M.alma)).astype(int)
end = time.time()
print("This step lasted {0}s".format(end-start))

print("Setting the permeability")
M.permeability[:] = 1
area = dx*dy

for i in range(len(M.coarse_volumes)):
    print("Assembly of coarse volume {0}".format(i))
    start = time.time()
    coef = lil_matrix((125, 125), dtype=np.float_)
    #for a in range(len(M.coarse_volumes)):
    for b in range(125):
        adjacencies = M.coarse_volumes[i].volumes.bridge_adjacencies(b, 2, 3)
        length = np.shape(adjacencies)
        for c in range(length[1]):
            id = np.array([adjacencies[0][c]],  dtype= np.int)
            coef[b,id] = equiv_perm(M.permeability[M.coarse_volumes[i].volumes.global_id[b]], M.permeability[M.coarse_volumes[i].volumes.global_id[id]])/centroid_dist(M.volumes.center[M.coarse_volumes[i].volumes.global_id[b]], M.volumes.center[M.coarse_volumes[i].volumes.global_id[id]])
        coef[b,b] = (-1)*coef[b].sum()
    end = time.time()
    print("This step lasted {0}s".format(end-start))

    print("Setting boundary conditions of coarse volume {0}".format(i))
    start = time.time()
    q = lil_matrix((125, 1), dtype=np.float_)
    coef[0:25] = 0
    q [0:25] = 500
    coef[100:125] = 0
    for r in range(25):
        coef[r,r] = 1
        coef[r+100,r+100] = 1
    end = time.time()
    print("This step lasted {0}s".format(end-start))

    print("Solving the problem of coarse volume {0}".format(i))
    start = time.time()
    coef = lil_matrix.tocsr(coef)
    q = lil_matrix.tocsr(q)
    P_coarse_volume = spsolve(coef,q)
    end = time.time()
    print("This step lasted {0}s".format(end-start))

    print("Storing results of coarse volume {0}".format(i))
    start = time.time()
    for a in range(125):
        M.coarse_volumes[i].pressure_coarse[a] = P_coarse_volume[a]
        #M.erro[i] = P[i]-P_analitico[i]
    end = time.time()
    print("This step lasted {0}s".format(end-start))

    for v in range(50):
        flux_rate = flux_rate + equiv_perm()*area*(M.coarse_volumes[])

print("Printing results")
M.core.print()
