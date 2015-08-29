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
import shutil
import copy
from datetime import datetime
from uuid import uuid4
from config import SCRATCH_DIR

class SlaterKosterXX:
    """ Class to manipulate the Hubbard U parameter
        in slater koster files of the type XX.spl - i.e.
        only ONE atom type.
    """

    def parse_slko_file(self, filename):

        f = open(filename, "r")
        lines = f.readlines()
        f.close()

        self.header = lines[0]
        self.footer = lines[2:]

        self.middle = lines[1]
        self.middle_left = "".join(lines[1].split(",")[:-2]) + ", "
        self.middle_right = "," + lines[1].split(",")[-1]

        self.Uhubb = float(self.middle.split()[4]) 

    def __init__(self, filename):

         self.parse_slko_file(filename)


    def write_slko_file(self, output_filename):

        output = copy.deepcopy(self.header)
        output += self.middle_left + "%11.6f  %11.6f  %11.6f" % \
            (self.Uhubb, self.Uhubb, self.Uhubb) + self.middle_right
        for foot in self.footer:
            output += foot

        f = open(output_filename, "w")
        f.write(output)
        f.close()

    def set_uhubb(self, uhubb_new):

        self.Uhubb = float(uhubb_new)



class ParameterSet:

    def _get_time_string(self):
        d = datetime.now()
        return str(d).replace(" ","_")

    def _generate_scr_dir(self):
        scr = SCRATCH_DIR + "/" + self._get_time_string() + "_" + str(uuid4())
        return scr

    def __init__(self, original_slko_dir, scr_dir=SCRATCH_DIR, zeta=None,
            Uhubb=None, Uhubb_derivatives=None):

        if zeta is None:
            self.zeta = 4.00
        else:
            self.zeta = zeta

        if Uhubb is None:
            # 3OB hubbard values
            self.Uhubbs = dict()
            self.Uhubbs["H"] = 0.4195 
            self.Uhubbs["C"] = 0.3647
            self.Uhubbs["N"] = 0.4309 
            self.Uhubbs["O"] = 0.4954
            self.Uhubbs["S"] = 0.3288
        else:
            self.Uhubb = Uhubb

        if Uhubb_derivatives is None:
            # 3OB hubbard values
            self.Uhubb_derivatives = dict()
            self.Uhubb_derivatives["H"] = -0.1857
            self.Uhubb_derivatives["C"] = -0.1492
            self.Uhubb_derivatives["N"] = -0.1535
            self.Uhubb_derivatives["O"] = -0.1575
            self.Uhubb_derivatives["S"] = -0.1100
        else:
            self.Uhubb_derivatives = Uhubb_derivatives

        self.angular_momentum = dict()
        self.angular_momentum["H"]  = "s"
        self.angular_momentum["C"]  = "p"
        self.angular_momentum["N"]  = "p"
        self.angular_momentum["O"]  = "p"
        self.angular_momentum["S"]  = "d"

        self.slko_filenames = dict()
        self.slko_filenames["H"] = "hh.spl"
        self.slko_filenames["C"] = "cc.spl"
        self.slko_filenames["N"] = "nn.spl"
        self.slko_filenames["O"] = "oo.spl"
        self.slko_filenames["S"] = "ss.spl"

        self.slko_tables = dict()
        self.slko_tables["H"] = SlaterKosterXX(original_slko_dir + "/hh.spl")
        self.slko_tables["C"] = SlaterKosterXX(original_slko_dir + "/cc.spl")
        self.slko_tables["N"] = SlaterKosterXX(original_slko_dir + "/nn.spl")
        self.slko_tables["O"] = SlaterKosterXX(original_slko_dir + "/oo.spl")
        self.slko_tables["S"] = SlaterKosterXX(original_slko_dir + "/ss.spl")

        self.slko_dir = self.make_slko_scr(original_slko_dir) + "/"


    def set_hubbard(self, atom_type, value):
        self.Uhubbs[atom_type] = value
        self.slko_tables[atom_type].set_uhubb(value)
        self.slko_tables[atom_type].write_slko_file(self.slko_dir \
                + self.slko_filenames[atom_type])

        # write new slater-koster file

    def set_hubbard_derivative(self, atom_type, value):
        self.Uhubb_derivatives[atom_type] = value

    def set_zeta(self, value):
        self.zeta = value

    def make_slko_scr(self, initial_slko_dir):

        slko_dir = self._generate_scr_dir()
        
        shutil.copytree(initial_slko_dir, slko_dir)
        # os.system("ls " + slko_dir)

        return slko_dir
