import os
from copy import deepcopy
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

def typeval_to_string(typeval):

    inputval = deepcopy(typeval)

    output = ""

    for key in ["S", "O", "N", "C", "H"]:

        n = inputval // TYPEVALS[key]

        inputval -= n * TYPEVALS[key]
        if n > 0:
            output += "%1s%1i " %(key, n)

    return "%-12s" % output



def get_typeval(obatom):

    name = obatom.GetType()[0]

    return TYPEVALS[name]

if __name__ == "__main__":

    np.set_printoptions(formatter={'float': '{: 0.3f}'.format}, linewidth=1000000)
    gaussian_mulliken = load_pickle("charges_gaussian.pickle")
    # charmm_mulliken = load_pickle("charges_3ob_npa.pickle")
    charmm_mulliken = load_pickle("charges_test.pickle")

    path = "xyz_sorted/"

    stats = dict()

    for atom in ATOMS:
        stats[atom] = dict()

    listing = os.listdir(path)
    for filename in sorted(listing):
        if filename.endswith(".xyz"):

            logfile = filename.replace(".xyz", ".log")
            dftb_mulliken = charmm_mulliken[filename]
            pbe_mulliken = gaussian_mulliken[logfile]

            qdiff = np.array(dftb_mulliken) - np.array(pbe_mulliken)

            max_qdiff = max(qdiff.min(), qdiff.max(), key=abs)
            # print
            # print "%-30s  %7.4f" % (filename, max_qdiff), qdiff
            # print "%39s" % "DFTB3/3OB", np.array(dftb_mulliken)
            # print "%39s" % "PBE/aug-cc-pVTZ", np.array(pbe_mulliken)

            mol = pybel.readfile("xyz", path + filename).next()

            for i, atom in enumerate(mol):
                type_int = 0

                # print "%-6s bonds to: " % (atom.OBAtom.GetType()),

                bonds = ""
                for obatom in openbabel.OBAtomAtomIter(atom.OBAtom):
                    bonds += "%-6s" % obatom.GetType()
                    type_int += get_typeval(obatom)



                while len(bonds) < 24:
                    bonds += "-     "

                # print "%-24s" % bonds,
                # print "%7.3f   ID: %05i" % (qdiff[i], type_int)

                name = atom.OBAtom.GetType()[0]

                if type_int not in stats[name].keys():
                    stats[name][type_int] = []

                stats[name][type_int].append(qdiff[i])

    for atom in ATOMS:

        all_values = []

        for bond in sorted(stats[atom]):

            values = np.array(stats[atom][bond])

            typeval_to_string(bond)
            print "%2s   --   %12s   %7.4f  %7.4f  %3i" % (atom, typeval_to_string(bond), np.mean(values), np.std(values),
                    len(values))

            all_values += stats[atom][bond]

        if len(all_values) < 1:
            continue

        rmsd = np.sqrt(np.mean(np.square(np.array(all_values))))
        mean = np.mean(np.array(all_values))

        print "%s:  RMSD = %4.2f   mean = %4.2f" % (atom, rmsd, mean)




