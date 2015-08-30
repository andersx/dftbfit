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

import numpy as np



class Sampler:


    move_names = dict()
    move_names[0]  = "Zeta"
    move_names[1]  = "Uh"
    move_names[2]  = "Uhd"
    move_names[3]  = "Uc"
    move_names[4]  = "Ucd"
    move_names[5]  = "Un"
    move_names[6]  = "Und"
    move_names[7]  = "Uo"
    move_names[8]  = "Uod"
    move_names[9]  = "Us"
    move_names[10] = "Usd"


    def __init__(self, par):

        self.par = par

        # Some initial values for these parameters
        self.move = 0
        self.change = 1.0

        self.move_old = 0
        self.change_old = 1.0

 
    def _change_parameters(self, move, change):

        if move == 0:
            zeta = self.par.zeta*change
            self.par.set_zeta(zeta)
        if move == 1:
            hubbard = self.par.Uhubbs["H"]*change
            self.par.set_hubbard("H", hubbard)
        if move == 2:
            hubbard_derivative = self.par.Uhubb_derivatives["H"]*change
            self.par.set_hubbard_derivative("H", hubbard_derivative)
        if move == 3:
            hubbard = self.par.Uhubbs["C"]*change
            self.par.set_hubbard("C", hubbard)
        if move == 4:
            hubbard_derivative = self.par.Uhubb_derivatives["C"]*change
            self.par.set_hubbard_derivative("C", hubbard_derivative)
        if move == 5:
            hubbard = self.par.Uhubbs["N"]*change
            self.par.set_hubbard("N", hubbard)
        if move == 6:
            hubbard_derivative = self.par.Uhubb_derivatives["N"]*change
            self.par.set_hubbard_derivative("N", hubbard_derivative)
        if move == 7:
            hubbard = self.par.Uhubbs["O"]*change
            self.par.set_hubbard("O", hubbard)
        if move == 8:
            hubbard_derivative = self.par.Uhubb_derivatives["O"]*change
            self.par.set_hubbard_derivative("O", hubbard_derivative)
        if move == 9:
            hubbard = self.par.Uhubbs["S"]*change
            self.par.set_hubbard("S", hubbard)
        if move == 10:
            hubbard_derivative = self.par.Uhubb_derivatives["S"]*change
            self.par.set_hubbard_derivative("S", hubbard_derivative)

    def get_parameters(self):

        return self.par


    def sample(self):

        self.move = np.random.choice(range(11))
        # self.move = 0 #  np.random.choice(range(10))
        self.change = np.random.uniform(0.95, 1.05)
        self._change_parameters(self.move, self.change)

        return self.move_names[self.move], self.change

    def reject(self):

        self._change_parameters(self.move, 1.0/self.change)

    def iprint(self):


        print "-- Current parameters --------------------"
        print "Zeta:          %7.4f" % self.par.zeta

        for a in ["H", "C", "N", "O", "S"]:

            print "Hubbard U:     %7.4f" % self.par.Uhubbs[a]

        for a in ["H", "C", "N", "O", "S"]:
            print "Hubbard Ude:   %7.4f" % self.par.Uhubb_derivatives[a]
        
        print "------------------------------------------"





