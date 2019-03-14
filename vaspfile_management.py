#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 13:25:11 2019

@author: mitchell koerner
@email:  mfkoerner@gmail.com
"""




from pymatgen.io.vasp import Incar
from pymatgen.io.vasp.outputs import Procar, Poscar
from pymatgen.core import periodic_table as pt
from shutil import copy
from os.path import join as j


with open('magset_z.txt') as f:
    MAGSET_Z = {int(i.rstrip('\n')) for i in f.readlines()}
with open('magset.txt') as f:
    MAGSET = {i.rstrip('\n') for i in f.readlines()}


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

    def unset_relaxation(self):
        """ Removes tags associated with relaxation runs """
        self.remove_tag('IBRION')
        self.remove_tag('EDIFFG')
        self.remove_tag('NSW')

    def set_static(self):
        """
        Remember to copy CONTCAR to POSCAR
        """
        self.unset_relaxation()
        self['NSW'] = 0
        self['SIGMA'] = 0.002
        self['LORBIT'] = 12

    def set_band(self):
        """
        Should come from a static run (so NSW should already be zero)
        Need to create correct KPOINTS
        Need to copy over CHGCAR from static run
        """
        self.unset_relaxation()
        self['ISMEAR'] = 0
        self['LORBIT'] = 12
        self['ICHARG'] = 11
        self['SIGMA'] = 0.002

    def set_magnetic(self, staticpath = '../static'):
        """
        should come from a static run that has both CHGCAR and wavecar
        Can be used in combination with set_band to create a band structure run
        """
        self.unset_relaxation()
        self['ISPIN'] = 2   # spin polarized calculation
        self['ICHARG'] = 1  # use CHGCAR from previous
        self['MAGMOM'] = self.get_POSCAR_mag_init(rundir = staticpath)


    def set_wavecar(self, value = True):
        """
        Turns on wavecar generation by default
        """
        self['LWAVE'] = value

    def increase_NBANDS(self, staticpath = '../static', factor = 2):
        """
        multiplies NBANDS from staticpath/PROCAR by factor
        """
        ogbands = self.get_NBANDS(staticpath)
        newbands = factor * ogbands
        self['NBANDS'] = newbands

    def get_NBANDS(self, rundir='.'):
        filepath = j(rundir, 'PROCAR')
        pro = Procar(filepath)
        return pro.nbands

    def get_POSCAR_atoms(self, rundir = '.'):
        filepath = j(rundir, 'POSCAR')
        pos = Poscar.from_file(filepath)
        atoms = repeated_list(pos.site_symbols, pos.natoms)
        return atoms

    def get_POSCAR_mag_init(self, rundir = '.', default_value = 3.0):
        atoms = self.get_POSCAR_atoms(rundir = rundir)
        is_magnetic = [i in MAGSET for i in atoms]
        mag_init = [i*default_value for i in is_magnetic]
        return mag_init

    def is_compound_magnetic(self, rundir = '.'):
        mag_init = self.get_POSCAR_mag_init(rundir = rundir)
        any_magnetic = sum(mag_init) > 0
        return any_magnetic

def copy_inputs(olddir, newdir):
    """
    copies POSCAR, POTCAR, INCAR, KPOINTS from olddir to newdir
    """
    copy(j(olddir, 'POSCAR'), j(newdir, 'POSCAR'))
    copy(j(olddir, 'POTCAR'), j(newdir, 'POTCAR'))
    copy(j(olddir, 'INCAR'), j(newdir, 'INCAR'))
    copy(j(olddir, 'KPOINTS'), j(newdir, 'KPOINTS'))

def ztostr(z):
    return pt.get_el_sp(z).name

def get_name(matid):
    Xz = matid[1:3]
    Az = matid[3:5]
    Bz = matid[5:7]
    X = ztostr(Xz)
    A = ztostr(Az)
    B = ztostr(Bz)
    return(r'{}$_3${}{}'.format(X, A, B))

def repeated_list(names, nums):
    out = []
    for i, num in enumerate(nums):
        for j in range(num):
            out.append(names[i])
    return out



