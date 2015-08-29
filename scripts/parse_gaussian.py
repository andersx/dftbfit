import os
import cPickle
from sys import argv

ATOMS = ["H", "C", "N", "O", "S", "P"]

def read_charges(filename):
    """ Reads NPA charges from an outputfile containing NPA/NBO output
    """

    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    start_line = 0
    end_line = 0

    for i, line in enumerate(lines):

        if "Atom  No    Charge         Core      Valence    Rydberg      Total" in line:
            start_line = i+2

        if "* Total *" in line:
            end_line = i - 1

    data = lines[start_line:end_line]

    charges = []
    types = []

    for atom in ATOMS:
        for entry in data: 
            if atom in entry:

                tokens = entry.split()
                charge = float(tokens[2]) 

                charges.append(charge)
                types.append(tokens[0])

    return charges, types



if __name__ == "__main__":

    path = "gaussian/"

    npa_charges = dict()
    npa_types = dict()

    listing = os.listdir(path)
    for filename in sorted(listing):
        if filename.endswith(".log"):

            charges, types = read_charges(path + filename)

            npa_charges[filename] = charges
            npa_types[filename] = types

            # for i in range(len(charges)):
            #    print types[i], i+1, charges[i]

    listing = os.listdir(path)
    for filename in sorted(listing):
        if filename.endswith(".log"):
            print filename, npa_charges[filename]


    f = open("charges_gaussian.pickle","wb")
    cPickle.dump(npa_charges, f, protocol=2)
    f.close()

    # f = open("types.pickle","wb")
    # cPickle.dump(npa_types, f, protocol=2)
    # f.close()

