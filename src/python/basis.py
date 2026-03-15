import numpy as np

class BasisFunction:

    def __init__(self, center, exponents, coefficients):

        # center of the orbital (atom position)
        self.center = np.array(center)

        # primitive Gaussian exponents
        self.exponents = exponents

        # contraction coefficients
        self.coefficients = coefficients

        # number of primitives
        self.nprimitive = len(exponents)
