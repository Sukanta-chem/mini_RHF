from basis import BasisFunction

# STO-3G parameters for H and He

basis_data = {
    "H": {"exponents": [3.42525091, 0.62391373, 0.16885540], "coefficients": [0.15432897, 0.53532814, 0.44463454] },
    "He": { "exponents": [6.36242139, 1.15892300, 0.31364979], "coefficients": [0.15432897, 0.53532814, 0.44463454] }
}

def build_basis(molecule):

    basis_functions = []
    for i, atom in enumerate(molecule.symbols):

        data = basis_data[atom]
        bf = BasisFunction( center = molecule.coords[i], exponents = data["exponents"], coefficients = data["coefficients"] )

        basis_functions.append(bf)

    return basis_functions
