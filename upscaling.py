# UPSCALLING OF STRUCTURED MESHES
import numpy as np
import time
import pdb
import xlsxwriter
from math import pi
from pymoab import rng, types
from tpfa.boundary_conditions import BoundaryConditions
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.linalg import spsolve
from preprocessor import M
import os

print("Initializating mesh")
os.system('python preprocessor.py')
dx, dy, dz = 1, 1, 1
nx, ny, nz = 25, 25,25
cx, cy, cz = 5, 5, 5
rx, ry, rz = 5, 5, 5
num_elements = nx*ny*nz
num_elements_coarse = rx*ry*rz

def equiv_perm(k1, k2):
    return (2*k1*k2)/(k1 + k2)

def centroid_dist(c1, c2):
    return ((c1-c2)**2).sum()

print("Setting the permeability")
M.permeability[:] = 1
area = dx*dy

for i in range(len(M.coarse_volumes)):
    print("Assembly of coarse volume {0}".format(i))
    start = time.time()
    adj = M.coarse_volumes[i].volumes.bridge_adjacencies(M.coarse_volumes[i].volumes.all, 2, 3) # IDs locais
    perm = M.permeability[M.coarse_volumes[i].volumes.global_id[M.coarse_volumes[i].volumes.all]]
    center = M.coarse_volumes[i].volumes.center[M.coarse_volumes[i].volumes.all]
    coef = lil_matrix((num_elements_coarse, num_elements_coarse), dtype=np.float_)
    for b in range(num_elements_coarse):
        adjacencies = adj[b] # Array de IDs locais
        for c in range(len(adjacencies)):
            id = np.array( [adjacencies[c]],  dtype= np.int)
            coef[b,id] = equiv_perm(perm[b], perm[id])/1
        coef[b,b] = (-1)*coef[b].sum()
    end = time.time()
    print("This step lasted {0}s".format(end-start))

    print("Setting boundary conditions of coarse volume {0}".format(i))
    start = time.time()
    q = lil_matrix((num_elements_coarse, 1), dtype=np.float_)
    coef[0:rx*ry] = 0
    q [0:rx*ry] = 500
    coef[(num_elements_coarse)-(rx*ry):num_elements_coarse] = 0
    for r in range(rx*ry):
        coef[r,r] = 1
        coef[r+(num_elements_coarse)-(rx*ry),r+(num_elements_coarse)-(rx*ry)] = 1
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
    for v in range(rx*ry):
        flow_rate =  + equiv_perm(perm[v], perm[v+rx*ry])*area*(M.coarse_volumes[i].pressure_coarse[v]-M.coarse_volumes[i].pressure_coarse[v+rx*ry])
        total_flow = total_flow + flow_rate

    permeability_coarse = total_flow/((area*rx*ry)*(M.coarse_volumes[i].pressure_coarse[v]-M.coarse_volumes[i].pressure_coarse[v+rx*ry]))
    print(permeability_coarse)

print("Assembly of upscaling")
start = time.time()
coef = lil_matrix((len(M.coarse_volumes), len(M.coarse_volumes)), dtype=np.float_)

#for i in range(len(M.coarse_volumes)):
#    coarse_centroid[i] = (sum(M.coarse_volumes[i].volumes.center)/125)

for i in range(len(M.coarse_volumes)):
    #M.coarse_volumes[i].permeability_coarse[:] = permeability_coarse
    #perm = M.coarse_volumes[i].permeability_coarse[:]
    adj = M.coarse_volumes[i].faces.coarse_neighbors
    for j in range(len(adj)):
        id = np.array(adj[j],  dtype= np.int)
        coef[i,id] = equiv_perm(1, 1)/25
    coef[i,i] = (-1)*coef[i].sum()
end = time.time()
print("This step lasted {0}s".format(end-start))

print("Setting boundary conditions of coarse mesh")
start = time.time()
q = lil_matrix((len(M.coarse_volumes), 1), dtype=np.float_)
coef[0:25] = 0
q [0:25] = 500
coef[100:125] = 0
for r in range(25):
    coef[r,r] = 1
    coef[r+100,r+100] = 1
end = time.time()
print("This step lasted {0}s".format(end-start))

print("Solving the problem")
start = time.time()
coef = lil_matrix.tocsr(coef)
q = lil_matrix.tocsr(q)
P = spsolve(coef,q)
end = time.time()
print("This step lasted {0}s".format(end-start))


for cvolume,index in zip(M.coarse_volumes,range(len(P))):
    M.pressure[cvolume.volumes.global_id[cvolume.volumes.all]] = P[index]


# for cvolume,index in zip(M.coarse_volumes,range(len(P))):
#     cvolume.pressure_coarse[:] = P[index]

#
# for i in range(len(M.coarse_volumes)):
#     M.pressure[M.coarse_volumes[i].volumes.global_idM.coarse_volumes[i].volumes.all] = P[i]

print("Printing results")
M.core.print()
