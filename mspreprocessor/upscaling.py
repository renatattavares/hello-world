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
nx, ny, nz = 20, 20, 20
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
    adj = M.coarse_volumes[i].volumes.bridge_adjacencies(M.coarse_volumes[i].volumes.all, 2, 3) # IDs locais
    perm = M.permeability[M.coarse_volumes[i].volumes.global_id[M.coarse_volumes[i].volumes.all]]
    center = M.coarse_volumes[i].volumes.center[M.coarse_volumes[i].volumes.all]
    coef = lil_matrix((125, 125), dtype=np.float_)
    for b in range(125):
        adjacencies = adj[b] # Array de IDs locais
        for c in range(len(adjacencies)):
            id = np.array( [adjacencies[c]],  dtype= np.int)
            coef[b,id] = equiv_perm(perm[b], perm[id])/centroid_dist(center[b], center[id])
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
    M.coarse_volumes[i].pressure_coarse[:] = P_coarse_volume
    end = time.time()
    print("This step lasted {0}s".format(end-start))

    total_flow = 0.0
    flow_rate = 0.0
    for v in range(25):
        flow_rate =  + equiv_perm(perm[v], perm[v+25])*area*(M.coarse_volumes[i].pressure_coarse[v]-M.coarse_volumes[i].pressure_coarse[v+25])
        total_flow = total_flow + flow_rate

    permeability_coarse = total_flow/((area*25)*(M.coarse_volumes[i].pressure_coarse[v]-M.coarse_volumes[i].pressure_coarse[v+25]))
    print(permeability_coarse)

print("Assembly of upscaling")
start = time.time()
center = M.volumes.center[M.volumes.all]
coef = lil_matrix((64, 64), dtype=np.float_)

for i in range(len(M.coarse_volumes)):
    #M.coarse_volumes[i].permeability_coarse[:] = permeability_coarse
    #perm = M.coarse_volumes[i].permeability_coarse[:]
    adj = M.coarse_volumes[i].faces.coarse_neighbors
    coarse_centroid = sum(M.coarse_volumes[i].volumes.center)/125
    for j in range(len(adj)):
        id = np.array(adj[j],  dtype= np.int)
        coarse_centroid_adj = sum(M.coarse_volumes[id].volumes.center)/125
        coef[i,id] = equiv_perm(1, 1)/centroid_dist(coarse_centroid, coarse_centroid_adj)
    coef[i,i] = (-1)*coef[i].sum()
end = time.time()
print("This step lasted {0}s".format(end-start))

print("Setting boundary conditions of coarse mesh")
start = time.time()
q = lil_matrix((64, 1), dtype=np.float_)
coef[0:16] = 0
q [0:16] = 500
coef[48:64] = 0
for r in range(16):
    coef[r,r] = 1
    coef[r+48,r+48] = 1
end = time.time()
print("This step lasted {0}s".format(end-start))

print("Solving the problem")
start = time.time()
coef = lil_matrix.tocsr(bc.coef)
q = lil_matrix.tocsr(bc.q)
P = spsolve(coef,q)
end = time.time()
print("This step lasted {0}s".format(end-start))

for i in range(len(M.coarse_volumes)):
    M.coarse_volumes[i].pressure_coarse[:] = P[i]

print("Printing results")
M.core.print()
