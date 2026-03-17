import numpy as np

def orthogonalization_matrix(S):
    eigvals, eigvecs = np.linalg.eigh(S)
    S_inv_sqrt = eigvecs @ np.diag(eigvals**(-0.5)) @ eigvecs.T

    return S_inv_sqrt


def build_density(C, nocc):
    n = C.shape[0]
    P = np.zeros((n, n))
    
    for mu in range(n):
        for nu in range(n):
            for m in range(nocc):
                P[mu, nu] += 2 * C[mu, m] * C[nu, m]

    return P


def build_fock(H, ERI, P):
    n = H.shape[0]
    F = np.zeros((n, n))

    for mu in range(n):
        for nu in range(n):
            F[mu, nu] = H[mu, nu]

            for lam in range(n):
                for sig in range(n):
                    F[mu, nu] += P[lam, sig] * ( ERI[mu, nu, lam, sig] - 0.5 * ERI[mu, sig, lam, nu] )

    return F


def compute_energy(P, H, F):
    E = 0.0
    n = P.shape[0]
    
    for mu in range(n):
        for nu in range(n):
            E += 0.5 * P[mu, nu] * (H[mu, nu] + F[mu, nu])

    return E


def scf(H, S, ERI, nocc, Enuc):
    X = orthogonalization_matrix(S)
    n = H.shape[0]
    P = np.zeros((n, n))
    E_old = 0.0

    for iteration in range(50):
        F = build_fock(H, ERI, P)

        # transform Fock
        Fp = X.T @ F @ X

        # diagonalize
        eps, Cp = np.linalg.eigh(Fp)

        # back transform
        C = X @ Cp

        # new density
        P = build_density(C, nocc)

        # energy
        E_elec = compute_energy(P, H, F)
        E_total = E_elec + Enuc
        print(f"Iter {iteration+1}: Energy = {E_total:.8f}")

        if abs(E_total - E_old) < 1e-8:
            print("\nSCF Converged")
            return E_total

        E_old = E_total

    print("SCF did not converge")
    return E_total
