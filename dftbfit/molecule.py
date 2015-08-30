#!/usr/bin/env python2
#
# Copyright (c) 2015, Anders Steen Christensen
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# 
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


import os
import copy
from datetime import datetime
from uuid import uuid4

from config import SCRATCH_DIR
from config import DFTB_EXE

class Molecule:

    def _get_atom_types(self, lines):

        atoms = []

        for line in lines[2:]:

            tokens = line.split()
            if not len(tokens) == 4:
                break
            
            atom = tokens[0]

            if atom not in atoms:
                atoms.append(atom)

        self.atom_types = atoms
        return atoms

    def _generate_atoms_dictionary(self, lines):

        atoms = self._get_atom_types(lines)

        atom_dictionary = dict()
        
        for i, atom_type in enumerate(atoms):

            atom_dictionary[atom_type] = i + 1

        return atom_dictionary

    def _print_header(self, lines):

        atoms = self._get_atom_types(lines)

        header = ""
        header += lines[0].split()[0] + "  C\n"

        header +=  " "
        for atom in atoms:
            header += " " + atom

        return header

    def _read_geometry(self, xyz_file):
        xyz_file = open(xyz_file)
        lines = xyz_file.readlines()
        xyz_file.close()

        atom_dictionary = self._generate_atoms_dictionary(lines)

        geometry = "Geometry = GenFormat {\n"
        geometry += self._print_header(lines)

        for i, line in enumerate(lines[2:]):

            tokens = line.split()


            geometry += "\n%4i %3i %16.10f %16.10f %16.10f " % \
                (i + 1 , atom_dictionary[tokens[0]], 
                float(tokens[1]),  float(tokens[2]), float(tokens[3]))

        geometry += "\n}"

        return geometry

    def _get_time_string(self):
        d = datetime.now()
        return str(d).replace(" ","_")

    def _generate_scr_dir(self):
        scr = SCRATCH_DIR + "/" + self._get_time_string() + "_" + str(uuid4())
        try:
            os.mkdir(scr)
        except:
            pass
        return scr

    def _read_number_of_atoms(self, xyz_file):

        xyz_file = open(xyz_file)
        lines = list(xyz_file.readlines())
        xyz_file.close()
        
        atoms = int(lines[0])

        return atoms



    def __init__(self, xyz_file, charge=0.0, 
            guess_charges=None, ref_charges=None):

        self.geometry = self._read_geometry(xyz_file)
        self.scr_dir = self._generate_scr_dir()
        self.charge = charge
        self.n_atoms = self._read_number_of_atoms(xyz_file)

        if guess_charges is None:
            self.guess_charges = [0.0 for _ in range(self.n_atoms)]

    def __del__(self):
        # os.system("ls " + self.scr_dir)
        try:
            os.remove(self.scr_dir + "/dftb_in.hsd")
        except:
            pass
        try:
            os.rmdir(self.scr_dir)
        except:
            pass

    def _write_dftbin(self, par, guess_charges=None):

        if guess_charges is None:
            guess_charges = [0.0 for _ in range(self.n_atoms)]

        output = copy.deepcopy(self.geometry)


        output += """\n\nHamiltonian = DFTB {
    charge = """
        output += str(self.charge)
        output += """\n    SCC = Yes
    SlaterKosterFiles {"""

        for atom_type1 in self.atom_types:
            for atom_type2 in self.atom_types:

                output += "\n        " + atom_type1.upper() + "-" + atom_type2.upper() + \
                        " = \"" + par.slko_dir + atom_type1.lower() + atom_type2.lower() + ".spl\""

        output += """\n    }
    MaxAngularMomentum {"""

        for atom_type in self.atom_types:
           output += "\n        %s = \"%s\" " % ( atom_type, par.angular_momentum[atom_type])

        output +=  """\n    }
    Filling = Fermi {
        Temperature [Kelvin] = 0.0
    }
    SCCTolerance = 1.0E-05
    Mixer = Broyden {}
    ThirdOrderFull = Yes
    MaxSCCIterations = 50
    DampXH = Yes"""
        output += "\n    DampXHExponent = " + str(par.zeta)
        output += """\n   HubbardDerivs = {
"""
        for atom_type in self.atom_types:
            output += "\n        %s = %s" % (atom_type, par.Uhubb_derivatives[atom_type])

        output += """\n    }

"""
        output += """
    InitialCharges = {
        AllAtomCharges = {"""
        for guess_charge in guess_charges:
            output += "        %f\n" % guess_charge
        output += """
        }
    }
}

Options {
    WriteBandOut = No
    WriteDetailedOut = No
    }

ParserOptions {
    ParserVersion = 4
    WriteHSDInput = No
}"""

        os.chdir(self.scr_dir)
        f = open("dftb_in.hsd", "w")
        f.write(output)
        f.close()
        # print output

    def run_dftb(self, par, verbose=False, guess_charges=None):

        if guess_charges is None:
            guess_charges = [0.0 for _ in range(self.n_atoms)]

        os.chdir(self.scr_dir)

        self._write_dftbin(par, guess_charges=guess_charges)
        output = os.popen(DFTB_EXE)
        output = list(output)

        if verbose:
            for line in output:
                print line,

        energy, charges, scf = self.parse_dftb_output(output)

        energy *= 627.509 # Hartree => kcal/mol

        return energy, charges, scf

    def parse_dftb_output(self, output):


        iscc_start = 0
        iscc_end = 0

        energy = 0.0

        imul_start = 0
        imul_end = 0

        for i, line in enumerate(output):

            if "iSCC Total electronic   Diff electronic      SCC error" in line:
                iscc_start = i+1
            if "Charges saved for restart in charges.bin" in line:
                iscc_end = i
            if "Net atomic charges (e) (begin)" in line:
                imul_start = i+1
            if "Net atomic charges (e) (end)" in line:
                imul_end = i
            if "Total Energy:" in line:
                energy = float(line.split()[2])

        scf = iscc_end - iscc_start

        charges = []
        for line in output[imul_start:imul_end]:
            charges.append(float(line.split()[1]))

        return energy, charges, scf

if __name__ == "__main__":

    import sys
    from parameterset import parameterset

    xyz_file = sys.argv[1]
    par = parameterset("/home/andersx/projects/nbo_test/hubbard_fit/slko_3OB/")

    mol = molecule(xyz_file)
    energy, charges, scf = mol.run_dftb(par, verbose=True)

    print mol.scr_dir
    print energy, charges, scf
