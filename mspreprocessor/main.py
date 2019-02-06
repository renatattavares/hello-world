
import numpy as np
import time
import pdb
import geoUtil.geoTools as gtool
import xlsxwriter
from math import pi, sqrt
from pymoab import rng, types
from meshHandle.multiscaleMesh import FineScaleMeshMS as msh
from meshHandle.corePymoab import CoreMoab as core
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.linalg import spsolve

################################################################################
def equiv_perm(k1, k2):
    return (2*k1*k2)/(k1 + k2)

def centroid_dist(c1, c2):
    return ((c1-c2)**2).sum()

def pressao_prescrita(coef, num_elements, nx, ny):
    q = lil_matrix((num_elements, 1), dtype=np.float_)
    coef[0:nx*ny] = 0
    q [0:nx*ny] = 500
    coef[num_elements-(nx*ny):num_elements] = 0
    for r in range(nx*ny):
        coef[r,r] = 1
        coef[r+num_elements-(nx*ny),r+num_elements-(nx*ny)] = 1
    return coef, q

################################################################################
dx, dy, dz = 1, 1, 1

print("Initializating mesh")
M = msh("malha-teste.h5m", dim = 3)
vec = np.arange(len(M.alma)).astype(int)

#################### Informações de condições de contorno ######################
num_bc = 1 # 1 para pressão prescrita, 2 para vazão, 3 para mista
value = 500
linha = 0
coluna = 0
nx = 10
ny = 10
nz = 10
num_elements = nx*ny*nz


################################################################################
print("Setting the variables")
M.permeability[:] = 1

print("Assembly")# Montagem da matriz de coeficientes do sistema.
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
coef, q = pressao_prescrita(coef, num_elements, nx, ny)

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
coef = lil_matrix.tocsr(coef)
q = lil_matrix.tocsr(q)
P = spsolve(coef,q)

print("Storing results")
for i in range(num_elements):
    M.pressure[i] = P[i]
    #M.erro[i] = P[i]-P_analitico[i]

print("Printing results")
M.core.print()
