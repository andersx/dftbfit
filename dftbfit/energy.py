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

import cPickle
import numpy as np

class Energy:

    def _load_pickle(self, filename):
        f = open(filename,"rb")
        p = cPickle.load(f)
        f.close()

        return p

    def __init__(self, reference_pickle):

        self.reference_data = self._load_pickle(reference_pickle)

        # Hack to clean up pickle :(
        for name in sorted(self.reference_data):
            new_name = name.rstrip(".log").rstrip(".xyz")
            self.reference_data[new_name] = self.reference_data[name]
            del self.reference_data[name]


    def calc_loglik(self, data, sigma=1.0, beta=1.0):

        n = 0
        chi_square = 0.0

        for name in sorted(data):

            mu = np.array(self.reference_data[name])
            x = np.array(data[name])
            try:
                chi_square += np.sum(np.square(mu - x))
            except:
                print name, data[name]
                return 100000.0, 0

            n += int(len(x))

        loglik = (n+1) * np.log(sigma) + chi_square / (2*sigma*sigma)
            
        rmsd = np.sqrt(chi_square/n)

        return loglik, n, rmsd

                

if __name__ == "__main__":


    energy = Energy("/home/andersx/projects/nbo_test/dftbfit/pickles/charges_gaussian.pickle")

    print energy.calc_loglik("lal")
