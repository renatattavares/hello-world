import numpy as np
import time
import xlsxwriter
from gerador_malha import GeradorMalha as gm

####################### Informações de entrada da malha ########################

nx = 5 # Número de elementos na direção x
ny = 5 # Número de elementos na direção y
nz = 5 # Número de elementos na direção Z
dx, dy, dz= 1.0, 1.0, 1.0 # Tamanho dos elementos nas direções x e y
dim = 2
num_elements = nx*ny*nz
################################################################################

# Inicializando a malha
malha = gm(nx,ny,nz,dx,dy,dz,num_elements)
print("Mesh created")
