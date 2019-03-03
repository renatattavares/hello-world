from gerador_malha import GeradorMalha as gm

####################### Informações de entrada da malha ########################

nx = 80 # Número de elementos na direção x
ny = 80 # Número de elementos na direção y
nz = 80 # Número de elementos na direção Z
dx, dy, dz= 1.0, 1.0, 1.0 # Tamanho dos elementos nas direções x e y
dim = 2
num_elements = nx*ny*nz
################################################################################

# Inicializando a malha
malha = gm(nx,ny,nz,dx,dy,dz,num_elements)
print("Mesh created")
