from molecule import Molecule
from sto3g import build_basis
from integrals import build_overlap_matrix
from integrals import build_kinetic_matrix
from integrals import build_nuclear_matrix
from integrals import build_eri_tensor
from scf import scf
from logger import Logger

mol = Molecule("../../examples/h2.xyz")
log = Logger("RHF.log")

log.write("---------------------------------------")
log.write("Restricted Hartree-Fock")
log.write("---------------------------------------\n")
log.write("Atoms: " + str(mol.symbols))

log.write("\nCoordinates:")
for c in mol.coords:
    log.write("  " + "  ".join(f"{x:10.6f}" for x in c))

energy = mol.nuclear_repulsion()

log.write("\nNuclear repulsion energy:")
log.write(f"{energy:.8f}")

basis = build_basis(mol)
log.write(f"\nNumber of basis functions: {len(basis)}")

for i, bf in enumerate(basis):
    log.write(f"\nBasis function {i}")
    log.write(f"center: {bf.center}")
    log.write(f"exponents: {bf.exponents}")


S = build_overlap_matrix(basis)
log.write("\nOverlap matrix:")
for row in S:
    log.write("  " + "  ".join(f"{x:10.6f}" for x in row))


T = build_kinetic_matrix(basis)
log.write("\nKinetic energy matrix:")
for row in T:
    log.write("  " + "  ".join(f"{x:10.6f}" for x in row))


V = build_nuclear_matrix(basis, mol)
log.write("\nNuclear attraction matrix:")
for row in V:
    log.write("  " + "  ".join(f"{x:10.6f}" for x in row))


H = T + V
log.write("\nCore Hamiltonian:")
for row in H:
    log.write("  " + "  ".join(f"{x:10.6f}" for x in row))


ERI = build_eri_tensor(basis)
log.write("\nSome ERI values:")
log.write(f"(0 0 | 0 0): {ERI[0,0,0,0]:.6f}")
log.write(f"(0 0 | 1 1): {ERI[0,0,1,1]:.6f}")
log.write(f"(0 1 | 0 1): {ERI[0,1,0,1]:.6f}")


nelec = 2
nocc = nelec // 2
E_total = scf(H, S, ERI, nocc, mol.nuclear_repulsion(), log)
log.write("\nFinal RHF Energy:")
log.write(f"{E_total:.8f} Hartree")
log.close()
