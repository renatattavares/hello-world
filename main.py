import numpy as np
import time
import xlsxwriter
from tpfa_moab.gerador_malha import GeradorMalha as gm
from tpfa_moab.condicoes_contorno import BoundaryConditions as bc
from pymoab import types, rng, topo_util
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.linalg import spsolve

################################################################################
def get_centroid_coords(elem):
    global dx, dy, dz
    v = malha.mbcore.get_coords(elem)
    centroid_x = v[1] + (dx/2)
    centroid_y = v[0] + (dy/2)
    centroid_z = v[2] + (dz/2)
    return np.array([centroid_x, centroid_y, centroid_z])

def equiv_perm(k1, k2):
    return (2*k1*k2)/(k1 + k2)

def centroid_dist(c1, c2):
    return (c1[0] - c2[0])**2 + (c1[1] - c2[1])**2 + (c1[2] - c2[2])**2

def pressao_prescrita(coef, num_elements, nx, ny):
    q = lil_matrix((num_elements, 1), dtype=np.float_)
    coef[0:nx*ny] = 0
    q [0:nx*ny] = 500
    coef[num_elements-(nx*ny):num_elements] = 0
    for r in range(nx*ny):
        coef[r,r] = 1
        coef[r+num_elements-(nx*ny),r+num_elements-(nx*ny)] = 1
    return coef, q

####################### Informações de entrada da malha ########################

nx = 20 # Número de elementos na direção x
ny = 20 # Número de elementos na direção y
nz = 20 # Número de elementos na direção Z
dx, dy, dz= 1.0, 1.0, 1.0 # Tamanho dos elementos nas direções x e y
dim = 2
num_elements = nx*ny*nz

#################### Informações de condições de contorno ######################
num_bc = 1 # 1 para pressão prescrita, 2 para vazão, 3 para mista
value = 500
linha = 0
coluna = 0

print("Initializating mesh")
# Inicializando a malha
malha = gm(nx,ny,nz,dx,dy,dz,num_elements)
mtu = topo_util.MeshTopoUtil(malha.mbcore)

print("Creating tags")
start = time.time()
# Determinando as coordenadas do centroide de cada elemento e armazenando-as em tags. Uma tag é um valor associado a cada elemento. Aqui, cada elemento possui duas tags: uma que armazena o valor das coordenadas do centroide e outra que armazena a permeabilidade.
centroid_tag = malha.mbcore.tag_get_handle('centroid', 3, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, True)
permeability_tag = malha.mbcore.tag_get_handle('permeability', 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, True)
pressure_tag = malha.mbcore.tag_get_handle('Pressão', 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, True)
id_tag = malha.mbcore.tag_get_handle('id', 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, True)
error_tag = malha.mbcore.tag_get_handle('erro', 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, True)
end = time.time()
print("This step lasted {0}s".format(end-start))

print("Setting data")
start = time.time()
i=0
for e in malha.elem_handles:
    elem_vertex = malha.mbcore.get_connectivity(e) # Função que informa quais os vérticies que compõem a entidade
    centroid_coord = get_centroid_coords(elem_vertex)
    malha.mbcore.tag_set_data(centroid_tag, e, centroid_coord)
    malha.mbcore.tag_set_data(id_tag, e, np.array([i],dtype=np.float_))
    if i <= 10:
        malha.mbcore.tag_set_data(permeability_tag, e, np.array([1], dtype=np.float_))
    else:
        malha.mbcore.tag_set_data(permeability_tag, e, np.array([1], dtype=np.float_))
    i = i+1
end = time.time()
print("This step lasted {0}s".format(end-start))

print("Assembly")# Montagem da matriz de coeficientes do sistema.
start = time.time()
coef = lil_matrix((num_elements, num_elements), dtype=np.float_)

for i in range(num_elements):
    bridge_adjacencies = mtu.get_bridge_adjacencies(malha.elem_handles[i], dim, 3, True)
    for j in range(len(bridge_adjacencies)):
        e1_tags = malha.mbcore.tag_get_tags_on_entity(bridge_adjacencies[j])
        e2_tags = malha.mbcore.tag_get_tags_on_entity(malha.elem_handles[i])
        e1_centroid = malha.mbcore.tag_get_data(e1_tags[0], bridge_adjacencies[j], flat=True)
        e2_centroid = malha.mbcore.tag_get_data(e2_tags[0], malha.elem_handles[i], flat=True)
        e1_perm = malha.mbcore.tag_get_data(e1_tags[1], bridge_adjacencies[j], flat=True)
        e2_perm = malha.mbcore.tag_get_data(e2_tags[1], malha.elem_handles[i], flat=True)
        e1_id = malha.mbcore.tag_get_data(e1_tags[2], bridge_adjacencies[j], flat=True)
        coef[i,e1_id] = equiv_perm(e1_perm, e2_perm)/centroid_dist(e1_centroid, e2_centroid)
    coef[i,i] = (-1)*coef[i].sum()
end = time.time()
print("This step lasted {0}s".format(end-start))

print("Setting boundary conditions")
start = time.time()
coef, q = pressao_prescrita(coef, num_elements, nx, ny)
end = time.time()
print("This step lasted {0}s".format(end-start))

'''
workbook = xlsxwriter.Workbook('coef_certo.xlsx')
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
coef = lil_matrix.tocsr(coef)
q = lil_matrix.tocsr(q)
P = spsolve(coef,q)
end = time.time()
print("This step lasted {0}s".format(end-start))

print("Printing results")
start = time.time()
i=0
for e in malha.elem_handles:
    malha.mbcore.tag_set_data(pressure_tag, e, P[i])
    #malha.mbcore.tag_set_data(error_tag, e, np.sqrt((P[i]-P_analitico[i])**2))
    i = i+1
malha.write_files()
end = time.time()
print("This step lasted {0}s".format(end-start))
