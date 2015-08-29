import os
import sys
import cPickle
import numpy as np
import openbabel
import pybel

sys.path.append('/home/andersx/dev/charmm-dftb-py')

from sccdftb_api import run_charmm, ATOMS

if __name__ == "__main__":

    path = "xyz_sorted/"

    npa_charges = dict()

    listing = os.listdir(path)
    for filename in sorted(listing):
        if filename.endswith(".xyz"):

            charge = 0.0
            if "[+1]" in filename:
                charge = 1.0
            elif "[-1]" in filename:
                charge = -1

            energy, dipole, scf, dftb_mulliken = run_charmm(path + filename,
                    clean_up=True, charge=charge, npa_charges=True, mixer=3,
                    scf_tol=1e-12, cpe=False)

            print "%-30s  %3i" % (filename, scf), dftb_mulliken
            npa_charges[filename] = dftb_mulliken

    f = open("charges_test.pickle","wb")
    cPickle.dump(npa_charges, f, protocol=2)
    f.close()
