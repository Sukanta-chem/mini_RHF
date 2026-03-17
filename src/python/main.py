from molecule import Molecule
from sto3g import build_basis
from integrals import build_overlap_matrix
from integrals import build_kinetic_matrix
from integrals import build_nuclear_matrix
from integrals import build_eri_tensor

mol = Molecule("../../examples/h2.xyz")

print("Atoms:")
print(mol.symbols)

print("\nCoordinates:")
print(mol.coords)

energy = mol.nuclear_repulsion()

print("\nNuclear repulsion energy:")
print(energy)

basis = build_basis(mol)
print("\nNumber of basis functions:", len(basis))

for i, bf in enumerate(basis):

    print("\nBasis function", i)
    print("center:", bf.center)
    print("exponents:", bf.exponents)


S = build_overlap_matrix(basis)
print("\nOverlap matrix:")
print(S)


T = build_kinetic_matrix(basis)
print("\nKinetic energy matrix:")
print(T)


V = build_nuclear_matrix(basis, mol)
print("\nNuclear attraction matrix:")
print(V)


H = T + V
print("\nCore Hamiltonian:")
print(H)


ERI = build_eri_tensor(basis)
print("\nSome ERI values:")
print("(0 0 | 0 0):", ERI[0,0,0,0])
print("(0 0 | 1 1):", ERI[0,0,1,1])
print("(0 1 | 0 1):", ERI[0,1,0,1])
