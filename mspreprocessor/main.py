
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

def pressao_prescrita(coef):
    q = lil_matrix((num_elements, 1), dtype=np.float_)
    coef[0:25] = 0
    q [0:25] = 500
    coef[100:125] = 0
    for r in range(25):
        coef[r,r] = 1
        coef[r+100,r+100] = 1
    return coef, q

def get_centroid_coords(elem):
    global dx, dy, dz
    v = M.core.mb.get_coords(elem)
    centroid_x = v[1] + (dx/2)
    centroid_y = v[0] + (dy/2)
    centroid_z = v[2] + (dz/2)
    return np.array([centroid_x, centroid_y, centroid_z])

P_analitico = np.array([500,500,500,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 ,500 , 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 375, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
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
num_elements = 5*5*5


################################################################################

print("Creating tags")
start = time.time()
# Determinando as coordenadas do centroide de cada elemento e armazenando-as em tags. Uma tag é um valor associado a cada elemento. Aqui, cada elemento possui duas tags: uma que armazena o valor das coordenadas do centroide e outra que armazena a permeabilidade.
centroid_tag = M.core.mb.tag_get_handle('centroid', 3, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, True)
permeability_tag = M.core.mb.tag_get_handle('permeability', 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, True)
pressure_tag = M.core.mb.tag_get_handle('Pressão', 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, True)
id_tag = M.core.mb.tag_get_handle('id', 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, True)
error_tag = M.core.mb.tag_get_handle('erro', 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, True)
end = time.time()
print("This step lasted {0}s".format(end-start))

print("Setting data")
start = time.time()
i=0
for e in M.volumes.elements_handle:
    elem_vertex = M.core.mb.get_connectivity(e) # Função que informa quais os vérticies que compõem a entidade
    centroid_coord = get_centroid_coords(elem_vertex)
    M.core.mb.tag_set_data(centroid_tag, e, centroid_coord)
    M.core.mb.tag_set_data(id_tag, e, np.array([i],dtype=np.float_))
    if i <= 73:
        M.core.mb.tag_set_data(permeability_tag, e, np.array([10], dtype=np.float_))
    else:
        M.core.mb.tag_set_data(permeability_tag, e, np.array([1], dtype=np.float_))
    i = i+1
end = time.time()
print("This step lasted {0}s".format(end-start))

print("Assembly")# Montagem da matriz de coeficientes do sistema.
start = time.time()
coef = lil_matrix((num_elements, num_elements), dtype=np.float_)

for i in range(num_elements):
    adjacencies = M.core.mtu.get_bridge_adjacencies(M.volumes.elements_handle[i], 2, 3, True)
    for j in range(len(adjacencies)):
        e1_tags = M.core.mb.tag_get_tags_on_entity(adjacencies[j])
        e2_tags = M.core.mb.tag_get_tags_on_entity(M.volumes.elements_handle[i])
        e1_centroid = M.core.mb.tag_get_data(e1_tags[0], adjacencies[j], flat=True)
        e2_centroid = M.core.mb.tag_get_data(e2_tags[0], M.volumes.elements_handle[i], flat=True)
        e1_perm = M.core.mb.tag_get_data(e1_tags[1], adjacencies[j], flat=True)
        e2_perm = M.core.mb.tag_get_data(e2_tags[1], M.volumes.elements_handle[i], flat=True)
        e1_id = M.core.mb.tag_get_data(e1_tags[2], adjacencies[j], flat=True)
        coef[i,e1_id] = equiv_perm(e1_perm, e2_perm)/centroid_dist(e1_centroid, e2_centroid)
    coef[i,i] = (-1)*coef[i].sum()
end = time.time()
print("This step lasted {0}s".format(end-start))

#print("Setting the permeability")
#M.permeability[:] = np.array([1])

#print("Assembly")# Montagem da matriz de coeficientes do sistema.
#coef = lil_matrix((num_elements, num_elements), dtype=np.float_)
#for i in range(num_elements):
#    adjacencies = M.volumes.bridge_adjacencies(i, 2, 3)
#    for j in range(len(adjacencies)):
#        coef[i,M.volumes.global_id[j]] = equiv_perm(M.permeability[i], M.permeability[M.volumes.global_id[j]])/centroid_dist(M.volumes.center[i], M.volumes.center[M.volumes.global_id[j]])
#    coef[i,i] = (-1)*coef[i].sum()

print("Setting boundary conditions")
coef, q = pressao_prescrita(coef)

print("Solving the problem")
coef = lil_matrix.tocsr(coef)
q = lil_matrix.tocsr(q)
P = spsolve(coef,q)

print("Storing results")
pressure_tag = M.core.mb.tag_get_handle('Pressure', 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, True)
error_tag = M.core.mb.tag_get_handle('Erro', 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, True)
i=0
for e in M.volumes.elements_handle:
    M.core.mb.tag_set_data(pressure_tag, e, P[i])
    M.core.mb.tag_set_data(error_tag, e, (P[i]-P_analitico[i]))
    i = i+1

print("Printing results")
M.core.print()
