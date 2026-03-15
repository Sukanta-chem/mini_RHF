from molecule import Molecule

mol = Molecule("../../examples/h2.xyz")

print("Atoms:")
print(mol.symbols)

print("\nCoordinates:")
print(mol.coords)

energy = mol.nuclear_repulsion()

print("\nNuclear repulsion energy:")
print(energy)
