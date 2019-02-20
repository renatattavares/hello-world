import numpy as np
from scipy.sparse import csr_matrix, lil_matrix

class BoundaryConditions():
    def __init__(self, num_elements, nx, ny, coef):
        self.coef = coef
        self.num_elements = num_elements
        self.nx = nx
        self.ny = ny
        self.coef, self.q = self.pressao_prescrita()

    def pressao_prescrita(self):
        self.q = lil_matrix((self.num_elements, 1), dtype=np.float_)
        self.coef[0:self.nx*self.ny] = 0
        self.q [0:self.nx*self.ny] = 500
        self.coef[self.num_elements-(self.nx*self.ny):self.num_elements] = 0
        for r in range(self.nx*self.ny):
            self.coef[r,r] = 1
            self.coef[r+self.num_elements-(self.nx*self.ny),r+self.num_elements-(self.nx*self.ny)] = 1
        return self.coef, self.q
