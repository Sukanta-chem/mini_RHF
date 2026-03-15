import numpy as np

class Molecule:

    def __init__(self, xyz_file):

        self.symbols = []
        self.coords = []

        self.read_xyz(xyz_file)
        self.coords = np.array(self.coords)

    def read_xyz(self, xyz_file):

        with open(xyz_file) as f:
            lines = f.readlines()

        lines = lines[2:]

        for line in lines:
            parts = line.split()
            atom = parts[0]

            x = float(parts[1])
            y = float(parts[2])
            z = float(parts[3])

            self.symbols.append(atom)
            self.coords.append([x, y, z])

    def nuclear_repulsion(self):

        charge = {"H":1, "He":2}
        energy = 0.0
        n = len(self.symbols)

        for i in range(n):
            Zi = charge[self.symbols[i]]

            for j in range(i+1, n):
                Zj = charge[self.symbols[j]]

                r = np.linalg.norm(self.coords[i] - self.coords[j])
                energy += Zi * Zj / r

        return energy
