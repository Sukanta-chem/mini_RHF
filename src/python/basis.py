import numpy as np

class BasisFunction:

    def __init__(self, center, exponents, coefficients):

        self.center = np.array(center)
        self.exponents = exponents
        self.coefficients = coefficients
        self.nprimitive = len(exponents)
