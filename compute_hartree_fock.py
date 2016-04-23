from hartree_fock import HartreeFock
from molecular_geometry import GeomReader
import profile
import sys


def compute_hartree_fock(mol_file,bas_file):
    geom = GeomReader(mol_file)
    mol = geom.mol
    h = HartreeFock(mol,bas_file,1)
    h.initial_calculation()
    h.converge(h.F, h.H_core, h.G, h.P,h.C,h.S_arb,h.ERI,h.aof,h.mol)
    return

if __name__=='__main__':
    if (len(sys.argv) <3):
       print 'Usage %s  [molecule_file] [basis_set_file]'%sys.argv[0]
       exit(0)
    compute_hartree_fock(sys.argv[2],sys.argv[3])
