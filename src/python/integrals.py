import numpy as np
import math

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


def boys_function(T):
    if T < 1e-8:
        return 1.0
    else:
        return 0.5 * np.sqrt(np.pi / T) * math.erf(np.sqrt(T))


def primitive_nuclear(alpha, beta, A, B, C, Z):

    R2 = np.dot(A - B, A - B)
    p = alpha + beta
    P = (alpha * A + beta * B) / p
    RPC2 = np.dot(P - C, P - C)
    expo = np.exp(-alpha * beta / p * R2)
    
    T = p * RPC2
    F0 = boys_function(T)
    V = -Z * (2 * np.pi / p) * expo * F0

    return V


def nuclear_attraction(bf1, bf2, molecule):

    Vtotal = 0.0
    charges = {"H":1, "He":2}

    for atom_index, atom in enumerate(molecule.symbols):
        C = molecule.coords[atom_index]
        Z = charges[atom]

        for i in range(bf1.nprimitive):
            alpha = bf1.exponents[i]
            ci = bf1.coefficients[i]
            Ni = (2*alpha/np.pi)**0.75

            for j in range(bf2.nprimitive):
                beta = bf2.exponents[j]
                cj = bf2.coefficients[j]
                Nj = (2*beta/np.pi)**0.75

                Vij = primitive_nuclear(alpha, beta, bf1.center, bf2.center, C, Z)
                Vtotal += ci * cj * Ni * Nj * Vij

    return Vtotal


def build_nuclear_matrix(basis, molecule):

    n = len(basis)
    V = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            V[i, j] = nuclear_attraction(basis[i], basis[j], molecule)

    return V


def primitive_eri(alpha, beta, gamma, delta, A, B, C, D):

    RAB2 = np.dot(A - B, A - B)
    RCD2 = np.dot(C - D, C - D)
    
    p = alpha + beta
    q = gamma + delta

    P = (alpha * A + beta * B) / p
    Q = (gamma * C + delta * D) / q

    RPQ2 = np.dot(P - Q, P - Q)

    K_ab = np.exp(-alpha * beta / p * RAB2)
    K_cd = np.exp(-gamma * delta / q * RCD2)
    
    T = (p * q / (p + q)) * RPQ2

    F0 = boys_function(T)
    prefactor = 2 * (np.pi ** 2.5) / (p * q * np.sqrt(p + q))

    return prefactor * K_ab * K_cd * F0


def eri(bf1, bf2, bf3, bf4):
    value = 0.0

    for i in range(bf1.nprimitive):
        alpha = bf1.exponents[i]
        ci = bf1.coefficients[i]
        Ni = (2*alpha/np.pi)**0.75

        for j in range(bf2.nprimitive):
            beta = bf2.exponents[j]
            cj = bf2.coefficients[j]
            Nj = (2*beta/np.pi)**0.75

            for k in range(bf3.nprimitive):
                gamma = bf3.exponents[k]
                ck = bf3.coefficients[k]
                Nk = (2*gamma/np.pi)**0.75

                for l in range(bf4.nprimitive):
                    delta = bf4.exponents[l]
                    cl = bf4.coefficients[l]
                    Nl = (2*delta/np.pi)**0.75

                    val = primitive_eri(alpha, beta, gamma, delta, bf1.center, bf2.center, bf3.center, bf4.center)
                    value += ci * cj * ck * cl * Ni * Nj * Nk * Nl * val

    return value


def build_eri_tensor(basis):

    n = len(basis)
    eri_tensor = np.zeros((n, n, n, n))

    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    eri_tensor[i, j, k, l] = eri( basis[i], basis[j], basis[k], basis[l] )

    return eri_tensor


