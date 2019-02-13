#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 13:25:11 2019

@author: mitchell koerner
@email:  mfkoerner@gmail.com
"""




from pymatgen.io.vasp import Incar

class PrefIncar(Incar):
    """
    Modifies the pymatgen Incar to include set functions that
     are for specific types of runs in Seshadri Group
    """

    def load_file(self, filepath = 'INCAR'):
        self.update(self.from_file(filepath))

    def set_system(self, system):
        self['System'] = system

    def remove_tag(self, tag):
        junk = self.pop(tag, None)

    def set_basics(self):
        self['GGA'] = 'Pe'
        self['LREAL'] = False
        self['EDIFF'] = 1e-6
        self['Prec'] = 'High'
        self['LWAVE'] = False
        self['LASPH'] = True
        self['ISMEAR'] = -5
        self['LMAXMIX'] = 6

    def set_full_relax(self):
        self['NSW'] = 30
        self['EDIFFG'] = -1e-3
        self['IBRION'] = 3
        self.remove_tag('SIGMA')

    def set_static(self):
        self['NSW'] = 0
        self['SIGMA'] = 0.002
        self['LORBIT'] = 12

    def set_band(self):
        """
        Should come from a static run (so NSW should already be zero)
        Need to create correct KPOINTS
        Need to copy over CHGCAR from static run
        """
        self['ISMEAR'] = 0
        self['LORBIT'] = 12
        self['ICHARG'] = 11
        self['SIGMA'] = 0.002








