from molecule import Molecule
from sto3g import build_basis
from integrals import build_overlap_matrix
from integrals import build_kinetic_matrix

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
