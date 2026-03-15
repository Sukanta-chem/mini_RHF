import numpy as np

def primitive_overlap(alpha, beta, A, B):

    # distance between centers
    R = np.linalg.norm(A - B)
    R2 = R * R

    p = alpha + beta
    prefactor = (np.pi / p) ** 1.5
    expo = np.exp(-alpha * beta / p * R2)

    return prefactor * expo

def overlap(bf1, bf2):

    S = 0.0

    for i in range(bf1.nprimitive):
        alpha = bf1.exponents[i]
        ci = bf1.coefficients[i]
        Ni = (2*alpha/np.pi)**0.75

        for j in range(bf2.nprimitive):
            beta = bf2.exponents[j]
            cj = bf2.coefficients[j]
            Nj = (2*beta/np.pi)**0.75

            Sij = primitive_overlap(alpha, beta, bf1.center, bf2.center)
            S += ci * cj * Ni * Nj * Sij

    return S

def build_overlap_matrix(basis):

    n = len(basis)
    S = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            S[i, j] = overlap(basis[i], basis[j])

    return S
