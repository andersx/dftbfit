import os
import sys
import cPickle
import numpy as np
import openbabel
import pybel

sys.path.append('/home/andersx/dev/charmm-dftb-py')

from sccdftb_api import run_charmm, ATOMS

def load_pickle(filename):
    f = open(filename,"rb")
    p = cPickle.load(f)
    f.close()
    return(p)

TYPEVALS = dict()
TYPEVALS["H"] = 1
TYPEVALS["C"] = 10
TYPEVALS["N"] = 100
TYPEVALS["O"] = 1000
TYPEVALS["S"] = 10000


def get_typeval(obatom):

    name = obatom.GetType()[0]

    return TYPEVALS[name]

if __name__ == "__main__":

    np.set_printoptions(formatter={'float': '{: 0.3f}'.format}, linewidth=1000000)
    gaussian_mulliken = load_pickle("charges_gaussian.pickle")
    # NPA gaussian_mulliken = load_pickle("charges_3ob_npa.pickle")
    # charmm_mulliken = load_pickle("charges_3ob.pickle")
    charmm_mulliken = load_pickle("charges_test.pickle")

    path = "xyz_sorted/"

    listing = os.listdir(path)
    for filename in sorted(listing):
        if filename.endswith(".xyz"):

            logfile = filename.replace(".xyz", ".log")
            dftb_mulliken = charmm_mulliken[filename]
            pbe_mulliken = gaussian_mulliken[logfile]
            # NPA pbe_mulliken = gaussian_mulliken[filename]

            qdiff = np.array(dftb_mulliken) - np.array(pbe_mulliken)

            max_qdiff = max(qdiff.min(), qdiff.max(), key=abs)
            print
            print "%-30s  %7.4f" % (filename, max_qdiff), qdiff
            print "%39s" % "DFTB3/3OB", np.array(dftb_mulliken)
            print "%39s" % "PBE/aug-cc-pVTZ", np.array(pbe_mulliken)

            mol = pybel.readfile("xyz", path + filename).next()

            for i, atom in enumerate(mol):
                type_int = 0

                print "%-6s bonds to: " % (atom.OBAtom.GetType()),

                bonds = ""
                for obatom in openbabel.OBAtomAtomIter(atom.OBAtom):
                    bonds += "%-6s" % obatom.GetType()
                    type_int += get_typeval(obatom)

                while len(bonds) < 24:
                    bonds += "-     "

                print "%-24s" % bonds,
                print "%7.3f   ID: %05i" % (qdiff[i], type_int)



