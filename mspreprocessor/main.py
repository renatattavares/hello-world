
import numpy as np
import time
import pdb
import geoUtil.geoTools as gtool
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

P_analitico = np.array([500,500,500,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 , 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
################################################################################

print("Inicializando a malha a partir do pre processador")
M = msh("malha-teste.h5m", dim = 3)
vec = np.arange(len(M.alma)).astype(int)

#################### Informações de condições de contorno ######################
num_bc = 1 # 1 para pressão prescrita, 2 para vazão, 3 para mista
value = 500
linha = 0
coluna = 0
num_elements = 5*5*5

################################################################################

adjacencies = [M.core.mb.get_adjacencies(e, 2, True) for e in M.volumes.elements_handle]# Encontrando adjacências para preencher a matriz de conectividade

# Incialização da matriz de conectividade. (Neste caso, a conectividade é em relação aos elementos, ou seja, quais elementos são vizinhos.)
connectivity = np.zeros((num_elements, num_elements), dtype=np.bool_)

print("Gerando matriz das conectividades (booleana)")
# Para cada adjacência diferente, verifica-se se existem uma fronteira compartilhada. Caso positivo, os dois elementos são vizinhos e isto é indicado em connectivity.
i, j = 0, 0
for a in adjacencies:
    for b in adjacencies:
        if b != a:
            intersection = rng.intersect(a, b)
            if not intersection.empty():
                connectivity[i][j] = 1
                connectivity[j][i] = 1
        j += 1
    j = 0
    i += 1

print("Definindo permeabilidade do meio")
M.permeability[:] = np.array([1])

print("Solving the problem")

coef = lil_matrix((num_elements, num_elements), dtype=np.float_)
for i in range(num_elements):
    for j in range(num_elements):
        # Se dois elementos são vizinhos e não se trata do mesmo elemento, então
        # são recuperados os valores das tags e calculado o valor do coeficiente.
        if connectivity[i,j] == True and i != j:
            coef[i,j] = equiv_perm(M.permeability[i], M.permeability[j])/centroid_dist(M.volumes.center[i], M.volumes.center[j])
    coef[i,i] = (-1)*coef[i].sum()

def pressao_prescrita(coef):
    q = lil_matrix((num_elements, 1), dtype=np.float_)
    coef[0:25] = 0
    q [0:25] = 500
    coef[100:125] = 0
    for r in range(25):
        coef[r,r] = 1
        coef[r+100,r+100] = 1
    return coef, q

coef, q = pressao_prescrita(coef)

coef = lil_matrix.tocsr(coef)
q = lil_matrix.tocsr(q)

P = spsolve(coef,q)

pressure_tag = M.core.mb.tag_get_handle('Pressure', 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, True)
error_tag = M.core.mb.tag_get_handle('Erro', 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, True)
i=0
for e in M.volumes.elements_handle:
    M.core.mb.tag_set_data(pressure_tag, e, P[i])
    M.core.mb.tag_set_data(error_tag, e, (P[i]-P_analitico[i]))
    i = i+1

print("Printing results")
M.core.print()
