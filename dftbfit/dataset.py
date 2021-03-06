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
from molecule import Molecule

class DataSet:


    def read_molecules_from_directory(self, xyz_dir):

        mols = dict()

        file_extension = ".xyz"
        listing = os.listdir(xyz_dir)

        for filename in sorted(listing):
            if filename.endswith(file_extension):

                charge = 0.0
                if "[+1]" in filename:
                    charge = 1.0
                elif "[-1]" in filename:
                    charge = -1

                mol = Molecule(xyz_dir + "/" + filename, charge=charge)

                name = filename.rstrip(file_extension)

                mols[name] = mol

        return mols

    def __init__(self, xyz_dir):

        self.molecules = self.read_molecules_from_directory(xyz_dir)



    def run_all(self, par, verbose=False, guess_charges=None):


        energy = dict()
        charges = dict()
        scf = dict()

        for name in sorted(self.molecules):

            if guess_charges is None:
                energy[name], charges[name], scf[name] = \
                        self.molecules[name].run_dftb(par, verbose=verbose)

            else:
                energy[name], charges[name], scf[name] = \
                        self.molecules[name].run_dftb(par, verbose=verbose,
                                guess_charges=guess_charges[name])

        return energy, charges, scf
