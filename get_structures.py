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
import qmpy         #quantum materials database interface
from qmpy.analysis.symmetry import WyckoffSite
from qmpy.analysis.symmetry import Spacegroup
import numpy
import scipy.stats
import os
import sys
from django.db.models import F

from qmpy import Structure

sg221 = Spacegroup.get(221)
wsite_a = WyckoffSite.get('a', sg221)
wsite_b = WyckoffSite.get('b', sg221)
wsite_c = WyckoffSite.get('c', sg221)
wsite_d = WyckoffSite.get('d', sg221)

def upgs():
    import gs
    reload(gs)
    return

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

def mode(s):
    """
    I don't know if nan_policy is producing weird results
    for elements without an electronegativity value
    """
    return scipy.stats.mode([i.species for i in s.atoms], nan_policy = 'omit')[0][0].__str__()

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
invper_min = [i for i in nonoxF if is_inv(i, 'min')]
#inverse perovskite, maximum requirement
invper_max = [i for i in nonoxF if is_inv(i, 'max')]
#inverse perovskite, bond valence requirement
# invper_bv =  [i for i in nonoxF if getv(i, check_inverse = True)]

# get both types of inverse perovskites
type1 = [i for i in invper_max if wsite_c in [j.wyckoff for j in i.sites]]
type2 = [i for i in invper_max if wsite_d in [j.wyckoff for j in i.sites]]


print('update')