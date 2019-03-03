#Condições de contorno, operações geométricas

import numpy as np

class BoundaryConditions():
    def __init__(self, num_bc, value, linha, coluna, num_elements, coef):
        self.coef = coef
        self.q = np.zeros([num_elements, 1])
        self.num_bc = num_bc
        self.value = value
        self.linha = linha
        self.coluna = coluna
        if self.num_bc == 1:
            self.pressao_prescrita()
        elif self.num_bc == 2:
            self.vazao_prescrita()
        else:
            self. pressao_prescrita()
            self.vazao_prescrita()

    def pressao_prescrita(self):
        self.coef[0:24] = 0
        self.coef[99:124] = 0
        for r in range(25):
            self.coef[r][r] = 1
            self.coef[r+99][r+99] = 1
            self.q [r] = 500
        self.coef[124][123] = 0

    def vazao_prescrita(self):
        self.q [4] = 20
