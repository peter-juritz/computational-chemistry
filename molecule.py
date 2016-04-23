class Molecule:
    n_atoms = 0
    bn = 0
    atoms = []
    bonds = []
    def add_atom(self,a):
        self.atoms.append(a);
        self.n_atoms+=1
        return
    def add_bond(self,a,b):
        self.bonds.append([a,b])
        self.bn+=1

    def __init__(self):
        self.atoms = []
        return



