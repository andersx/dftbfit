import openbabel
import pybel


filename = "xyz_sorted/113_HO-CH2CH2-OH[0].xyz"

mol = pybel.readfile("xyz", filename).next()

for atom in mol:
    print atom.OBAtom.GetType() 

    for obatom in openbabel.OBAtomAtomIter(atom.OBAtom):
        print obatom.GetType(),
