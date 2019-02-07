
import numpy as np
import time
import pdb
import geoUtil.geoTools as gtool
import xlsxwriter
from math import pi, sqrt
from pymoab import rng, types
from meshHandle.multiscaleMesh import FineScaleMeshMS as msh
from meshHandle.corePymoab import CoreMoab as core
from tpfa.assembly import Assembly
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

print("Assembly")
start = time.time()
coef = lil_matrix((num_elements, num_elements), dtype=np.float_)
for i in range(num_elements):
    adjacencies = M.volumes.bridge_adjacencies(i, 2, 3)
    length = np.shape(adjacencies)
    for j in range(length[1]):
        id = np.array([adjacencies[0][j]],  dtype= np.int)
        coef[i,id] = equiv_perm(M.permeability[i], M.permeability[id])/centroid_dist(M.volumes.center[i], M.volumes.center[id])
    coef[i,i] = (-1)*coef[i].sum()
end = time.time()
print("This step lasted {0}s".format(end-start))

print("Setting boundary conditions")
start = time.time()
bc = BoundaryConditions(num_elements, nx, ny, coef)
end = time.time()
print("This step lasted {0}s".format(end-start))

'''
workbook = xlsxwriter.Workbook('coef_mspreprocessor.xlsx')
worksheet = workbook.add_worksheet()
matrix = lil_matrix.toarray(coef)

row = 0
col = 0

for row in range(125):
  for col in range(125):
    worksheet.write(row, col, matrix[row][col])

workbook.close()
'''

print("Solving the problem")
start = time.time()
coef = lil_matrix.tocsr(bc.coef)
q = lil_matrix.tocsr(bc.q)
P = spsolve(coef,q)
end = time.time()
print("This step lasted {0}s".format(end-start))

print("Storing results")
start = time.time()
for i in range(num_elements):
    M.pressure[i] = P[i]
    #M.erro[i] = P[i]-P_analitico[i]
end = time.time()
print("This step lasted {0}s".format(end-start))

area = dx*dy*25
coef = lil_matrix((50, 50), dtype=np.float_)
#for a in range(len(M.coarse_volumes)):
for b in range(25):
    adjacencies = M.coarse_volumes[0].volumes.bridge_adjacencies(b, 2, 3)
    length = np.shape(adjacencies)
    for c in range(length[1]):
        id = np.array([adjacencies[0][c]],  dtype= np.int)
        coef[b,id] = equiv_perm(M.permeability[M.coarse_volumes[0].volumes.global_id[b]], M.permeability[M.coarse_volumes[0].volumes.global_id[id]])/centroid_dist(M.volumes.center[M.coarse_volumes[0].volumes.global_id[b]], M.volumes.center[M.coarse_volumes[0].volumes.global_id[id]])
    coef[b,b] = (-1)*coef[b].sum()

coef = lil_matrix.tocsr(bc.coef)
dp = lil_matrix((1,50), dtype=np.float_)
dp [:] = 500
dp = lil_matrix.tocsr(dp)

q = np.dot(dp,coef)

print("Printing results")
M.core.print()
