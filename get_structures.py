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



class StructureStats():
    """Docstring here

    Bugs: 
    Does not update structure if you try to, just 
    create a new instance if you want a new structure
    """
    def __init__(self, structs):
        """ structs should be a list of oqmd structure objects """
        self.structs = structs
        self.n = len(self.structs)
        self.done_get_valences = False
        self.done_invper_valences = False
        self._to_pymatgen()
        self._elstats()
        self._get_is_inv()
    def _to_pymatgen(self):
        """ updates  """
        self.lattices = [i.cell for i in self.structs]
        self.atoms = [s.atoms for s in self.structs]
        self.species = [[j.species for j in s.atoms] for s in self.structs]
        self.coordinates = [[i.coord for i in s.atoms] for s in self.structs]
        self.pystructs = [pymatgen.Structure(self.lattices[i], self.species[i],
         self.coordinates[i]) for i in range(self.n)]
    def _get_pmg_valence(self, pystruct):
        """ Runs work for pymatgen valence code"""
        try:
            vs = BV.get_valences(pystruct)
        except:
            vs = None
        return(vs)
    def _electro(self, element):
        return pymatgen.Specie(element).X
    def get_valences(self, output = True):
        """
        returns list of valences for each compound 
        compounds that fail will have None instead
        """
        if not self.done_get_valences:
            self.valences = [self._get_pmg_valence(i) for i in self.pystructs]
            self.done_get_valences = True
        if output:
            return self.valences
    def invper_valences(self, output = True):
        """
        Returns a list of inverse perovskites by X site valence charge > 0
        """
        if not self.done_get_valences:
            self.get_valences(output = False)
        if not self.done_invper_valences:
            self.invper_valences = []
            for i, valence in enumerate(self.valences):
                if not valence is None:
                    if scipy.stats.mode(valence, nan_policy = 'omit')[0][0] > 0:
                        self.invper_valences.append(self.structs[i])
            self.done_invper_valences = True
        if output:
            return self.invper_valences
    def _elstats(self):
        """
        string_species is pymatgen Specie object
        electronegativities is lists of electronegativities
        emin is electronegativities min
        emax is electronegativities max
        emode is electronegativities mode
        """
        self.string_species = [[pymatgen.Specie(a.element.__str__()) for a in atoms] for atoms in self.atoms]
        self.electronegativities = [[i.X for i in string_species] for string_species in self.string_species]
        self.emin = numpy.array([min(i) for i in self.electronegativities])
        self.emax = numpy.array([max(i) for i in self.electronegativities])
        self.emode = numpy.array([scipy.stats.mode(i, nan_policy = 'omit')[0][0] for i in self.electronegativities])
    def _get_is_inv(self, method = 'max'):
        """
        Checks if s is an inverse perovskite by method of electronegativity
        if X-site is not highest electronegativity, returns True, otherwise false
        if method = min, checks that X-site is minimum electronegativity
        """
        if method == 'max':
            self.is_inv = numpy.not_equal(self.emax, self.emode)
        elif method == 'min':
            self.is_inv = numpy.equal(self.emin, self.emode)
        else:
            raise ValueError('method must be either "min" or "max"')
        self.inverse_perovskites = [self.structs[i] for i in range(self.n) if self.is_inv[i]]






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

SG221 = Spacegroup.get(221)
WSITE_A = WyckoffSite.get('a', SG221)
WSITE_B = WyckoffSite.get('b', SG221)
WSITE_C = WyckoffSite.get('c', SG221)
WSITE_D = WyckoffSite.get('d', SG221)

class BetterStructure():
    """
    Adding more features to a qmpy structure
    initialises from a qmpy structure
    """
    def __init__(self, s):
        self.structure = s
        self.wyckoffsites = [i.wyckoff for i in self.structure.sites]
    def get_element_by_site(self, site):
        """
        gets atom at lettered site
        for given structure
        """
        spacegroup = self.structure.spacegroup
        wsite = WyckoffSite.get(site, spacegroup)
        element = self.structure.atoms[self.wyckoffsites.where(wsite)].element
        return element

    def get_elements_ordered_by_wyckoff_site(self, sites):
        """
        returns list of elements ordered by list of strings
        called sites that represents the wyckoff sites
        """
        ordered_elements = [self.get_element_by_site(site) for site in sites]
        return ordered_elements

    # def label_by_sites(s):
    #     """
    #     assigns label and id based on sites
    #     label goes X3AB
    #     id goes XXAABB
    #     may bug out if not simple cubic inverse perovskite
    #     may instead error out
    #     """
    #     if WSITE_C in self.wyckoffsites:
    #         ordered_elements = get_elements_ordered_by_wyckoff_site(s, ['a','b','c'])
    #     elif WSITE_D in self.wyckoffsites:

    #     else:
    #         raise ValueError("need simple cubic inverse perovskite")


# other stuff

allstructperovtheory = Structure.objects.filter(entry__meta_data__value__contains='perovskite', label='input')
icsd = Structure.objects.filter(
    entry__meta_data__value='icsd', 
    label='input',
    entry__id = F('entry__duplicate_of__id')
    )

# Get all unique icsd entries with 5 atoms and 3 unique species in spacegroup 221
un53icsd221 = Structure.objects.filter(
    entry__meta_data__value='icsd', 
    label='input',
    natoms = 5,
    ntypes = 3,
    entry__id = F('entry__duplicate_of__id'),
    spacegroup = 221
    )
# Knock out all with 3 Oxygen
nonox = [ i for i in un53icsd221 if not 'O3' in i.__str__() ]
# Knock out all with 3 Fluorine
nonoxF = [ i for i in nonox if not 'F3' in i.__str__() ]


nonoxFstr = [i.__str__() for i in nonoxF]

#inverse perovskite, minimum requirement
# invper_min = [i for i in nonoxF if is_inv(i, 'min')]
#inverse perovskite, maximum requirement
invper_max = [i for i in nonoxF if is_inv(i, 'max')]
#inverse perovskite, bond valence requirement
# invper_bv =  [i for i in nonoxF if getv(i, check_inverse = True)]

# filter out F block electrons
invper_nof = [i for i in invper_max if len(PARTIAL_F.intersection(set(i.elements))) == 0]

final_list = invper_nof
# get both types of inverse perovskites (searches for existance of certain wyckoff sites)
# type 1 A: a1; B: b1; X: c3; This is Ram's preferred way
type1 = {i for i in final_list if WSITE_C in [j.wyckoff for j in i.sites]}
# type 2 A: b1; B: a1; X: d3
type2 = {i for i in final_list if WSITE_D in [j.wyckoff for j in i.sites]}
# broken because of bug in qmpy.materials.structures.py definition of translate line 1462
type2_transformed = {struct.recenter(struct[[i.wyckoff.symbol for i in struct.sites].index(u'b')],
 in_place = False, middle = False) for struct in type2}





