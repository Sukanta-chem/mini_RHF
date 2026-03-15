import numpy as np

class BasisFunction:

    def __init__(self, center, exponents, coefficients):

        self.center = np.array(center)
        self.exponents = exponents
        self.coefficients = coefficients
        self.nprimitive = len(exponents)
        self.normalize()


    def primitive_overlap(self, alpha, beta):

        R2 = 0.0  
        p = alpha + beta
        prefactor = (np.pi / p) ** 1.5
        return prefactor * np.exp(-alpha * beta / p * R2)


    def normalize(self):

        S = 0.0

        for i in range(self.nprimitive):

            alpha = self.exponents[i]
            ci = self.coefficients[i]
            Ni = (2*alpha/np.pi)**0.75

            for j in range(self.nprimitive):

                beta = self.exponents[j]
                cj = self.coefficients[j]
                Nj = (2*beta/np.pi)**0.75

                Sij = self.primitive_overlap(alpha, beta)
                S += ci * cj * Ni * Nj * Sij

        norm = 1 / np.sqrt(S)

        self.coefficients = [c * norm for c in self.coefficients]
