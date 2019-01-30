#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 18:05:05 2019

@author: mitchell koerner
@email:  mfkoerner@gmail.com
"""


import pymatgen
from pymatgen.analysis.bond_valence import BVAnalyzer
BV = BVAnalyzer()
import inspect
# use with lines = inspect.getsource(foo)
import qmpy         #quantum materials database interface
from qmpy.analysis.symmetry import WyckoffSite
from qmpy.analysis.symmetry import Spacegroup
from qmpy.materials.element import Element
from qmpy import Structure
import numpy
import scipy.stats
import os
import sys
from django.db.models import F

# relative paths
import os
dirname = os.path.dirname(__file__)
def relpath(path):
    return os.path.join(dirname, path)



### important variables ###
# space group and wyckoff sites for inverse perovskites
SG221 = Spacegroup.get(221)
WSITE_A = WyckoffSite.get('a', SG221)
WSITE_B = WyckoffSite.get('b', SG221)
WSITE_C = WyckoffSite.get('c', SG221)
WSITE_D = WyckoffSite.get('d', SG221)
# partial f electrons
with open(relpath('data/partial_f.txt')) as f:
    PARTIAL_F_strs = [i.rstrip('\n') for i in f.readlines()]
PARTIAL_F = {Element.get(i) for i in PARTIAL_F_strs}


def mode(s):
    """
    I don't know if nan_policy is producing weird results
    for elements without an electronegativity value
    """
    return scipy.stats.mode([i.species for i in s.atoms], nan_policy = 'omit')[0][0].__str__()




# Old functions
# Convert oqmd structure to pymatgen structure
def pystruct(s):
    lattice = s.cell
    species = [i.species for i in s.atoms]
    coordinates = [i.coord for i in s.atoms]
    return pymatgen.Structure(lattice, species, coordinates)

def getv(s, check_inverse = False):
    """
    check_inverse only returns true if mode(valences) > 0
    """
    pys = pystruct(s)
    try:
        vs = BV.get_valences(pys)
    except:
        vs = None
    if check_inverse:
        if vs is None: return False
        elif scipy.stats.mode(vs, nan_policy = 'omit')[0][0] < 0:
            return False
        else: return True
    else:
        return vs

def electro(element):
    return pymatgen.Specie(element).X

def elstats(s):
    """
    Returns stats of electronegativity of s and electronegativity
    [min, max, mode, [electronegativity]]
    """
    sspec = [pymatgen.Specie(i.element.__str__()) for i in s.atoms]
    electronegativities = [i.X for i in sspec]
    emin = min(electronegativities)
    emax = max(electronegativities)
    emode = scipy.stats.mode(electronegativities, nan_policy = 'omit')[0][0]
    return [emin, emax, emode, electronegativities]
def is_inv(s, method = 'max'):
    """
    Checks if s is an inverse perovskite by method of electronegativity
    if X-site is not highest electronegativity, returns True, otherwise false
    if method = min, checks that X-site is minimum electronegativity
    """
    es = elstats(s)
    if method == 'max':
        return es[2] != es[1]
    elif method == 'min':
        return es[2] == es[0]
    else:
        raise ValueError('method must be either "min" or "max"')


class InversePerovskiteBonuses():
    """
    Adding more features to a qmpy structure
    initialises from a qmpy structure
    current main purpose is to label_by_sites

    ### bonuses ###
    wyckoffsites:       list of wyckoff sites
    ordered_elements:   elements in order X A B
    label:              X3AB
    idstr:              periodic numbers of order XXAABB
    """
    def __init__(self, s):
        self.structure = s
        self.wyckoffsites = [i.wyckoff for i in self.structure.sites]
        self._set_ordered_elements()
        self._set_label_by_sites()


    def _get_element_by_site(self, site):
        """
        gets atom at lettered site
        for given structure
        """
        spacegroup = self.structure.spacegroup
        wsite = WyckoffSite.get(site, spacegroup)
        element = self.structure.atoms[self.wyckoffsites.index(wsite)].element
        return element

    def _get_elements_ordered_by_wyckoff_sites(self, sites):
        """
        returns list of elements ordered by list of strings
        called sites that represents the wyckoff sites
        """
        ordered_elements = [self._get_element_by_site(site) for site in sites]
        return ordered_elements

    def _set_ordered_elements(self):
        """
        sets self.ordered elements to be in order X A B
        """
        if WSITE_C in self.wyckoffsites:
            self.ordered_elements = self._get_elements_ordered_by_wyckoff_sites(['c','a','b'])
        elif WSITE_D in self.wyckoffsites:
            self.ordered_elements = self._get_elements_ordered_by_wyckoff_sites(['d','b','a'])
        else:
            raise ValueError("need simple cubic inverse perovskite")

    def _set_label_by_sites(self):
        """
        Returns: label, idstr
        assigns label and id based on sites
        label goes X3AB
        id goes XXAABB
        may bug out if not simple cubic inverse perovskite
        may instead error out
        """
        self.label = '{}3{}{}'.format(*[i.symbol for i in self.ordered_elements])
        self.idstr = '{:02d}{:02d}{:02d}'.format(
            *[i.z for i in self.ordered_elements]
            )




