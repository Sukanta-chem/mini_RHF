import numpy as np

def primitive_overlap(alpha, beta, A, B):

    R = np.linalg.norm(A - B)
    R2 = np.dot(A - B, A - B)

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

def primitive_kinetic(alpha, beta, A, B):

    R2 = np.dot(A - B, A - B)
    p = alpha + beta
    prefactor = (np.pi / p) ** 1.5
    expo = np.exp(-alpha * beta / p * R2)

    S = prefactor * expo

    reduced = alpha * beta / p
    T = reduced * (3 - 2 * reduced * R2) * S

    return T


def kinetic(bf1, bf2):

    Tval = 0.0

    for i in range(bf1.nprimitive):
        alpha = bf1.exponents[i]
        ci = bf1.coefficients[i]
        Ni = (2*alpha/np.pi)**0.75

        for j in range(bf2.nprimitive):
            beta = bf2.exponents[j]
            cj = bf2.coefficients[j]
            Nj = (2*beta/np.pi)**0.75

            Tij = primitive_kinetic(alpha, beta, bf1.center, bf2.center)
            Tval += ci * cj * Ni * Nj * Tij

    return Tval


def build_kinetic_matrix(basis):

    n = len(basis)
    T = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            T[i, j] = kinetic(basis[i], basis[j])

    return T
