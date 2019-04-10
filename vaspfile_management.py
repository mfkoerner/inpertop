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


MAGSET_Z = {21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 39, 40, 41, 42, 43, 44, 
45, 46, 47, 48, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,
 72, 73, 74, 75, 76, 77, 78, 79, 80, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 
 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112}

MAGSET = {'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Y', 
'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'La', 'Ce', 'Pr', 
'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 
'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Ac', 'Th', 'Pa', 'U', 
'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 
'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn'}


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
        self['LWAVE'] = True
        self['LASPH'] = True
        self['ISMEAR'] = -5
        self['LMAXMIX'] = 6

    def start_from_CHGCAR(self):
        """
        Now this is used for magnetic and SOC instead of forcing it
        """
        self['ICHARG'] = 1

    def set_full_relax(self):
        """
        Always run this as the last move since many of the other ones
         specifically unset relaxation
        """
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
        self['LWAVE'] = True

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
        self['LWAVE'] = True
        self.remove_tag('MAGMOM')   # this should come from a preconverged CHGCAR

    def set_magnetic(self, rundir = '.', use_CHGCAR = True):
        """
        should come from a static run that has both CHGCAR and wavecar
        Can be used in combination with set_band to create a band structure run
        Make sure that POSCAR is copied before running set_magnetic
        THIS ONLY WORKS FOR COLINNEAR VASP, for noncolinnear use set_magnetic_ncl
        """
        self.unset_relaxation()
        if use_CHGCAR:
            self.start_from_CHGCAR()
        self['ISPIN'] = 2  # split polarization
        self['MAGMOM'] = self.get_POSCAR_mag_init(rundir = rundir)

    def set_magnetic_ncl(self, rundir = '.', use_CHGCAR = True):
        """
        This creates a magnetic noncolinnear run
        Remember to run this using vasp.ncl
        start from static run with both CHGCAR and WAVECAR if use_CHGCAR is true
        """
        self.unset_relaxation()
        if use_CHGCAR:
            self.start_from_CHGCAR()
        self['LNONCOLLINEAR'] = True
        self['ISPIN'] = 2
        self['MAGMOM'] = self.get_3d_POSCAR_mag_init(rundir = rundir)



    def set_SOC(self, use_CHGCAR = True):
        """
        Remember to run this using vasp.ncl
        Start from static run with both CHGCAR and WAVECAR
        """
        self.unset_relaxation()
        if use_CHGCAR:
            self.start_from_CHGCAR()
        self['LSORBIT'] = True
        self['LNONCOLLINEAR'] = True


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
        """
        Returns array of POSCAR atoms by order in POSCAR
        ex: ['A', 'X', 'X', 'X', 'B']
        """
        filepath = j(rundir, 'POSCAR')
        pos = Poscar.from_file(filepath)
        atoms = repeated_list(pos.site_symbols, pos.natoms)
        return atoms

    def get_POSCAR_mag_init(self, rundir = '.', default_value = 3.0):
        """
        returns array of 3.0 and 0.0 for magnetic vs nonmagnetic atoms in POSCAR
        """
        atoms = self.get_POSCAR_atoms(rundir = rundir)
        is_magnetic = [i in MAGSET for i in atoms]
        mag_init = [i*default_value for i in is_magnetic]
        return mag_init

    def get_3d_POSCAR_mag_init(self, rundir = '.', default_value = 3):
        """
        Returns 3d version of mag_init for use in noncollinear runs
        """
        mag_init = self.get_POSCAR_mag_init(rundir = rundir,
            default_value = default_value)
        mag_3d_init = [0.0 for i in range(len(mag_init) * 3)]
        for i in range(len(mag_init)):
            mag_3d_init[3*i] = mag_init[i]
        mag_3d_list_strings = ['{:3.1f}'.format(i) for i in mag_3d_init]
        mag_3d_string = ' '.join(mag_3d_list_strings)
        print(mag_3d_string)
        return mag_3d_string

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



