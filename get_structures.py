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
import numpy
import os
import sys
from django.db.models import F

from qmpy import Structure

def upgs():
    import gs
    reload(gs)
    return

allstruct = Structure.objects.filter(entry__meta_data__value__contains='perovskite', label='input')
icsd = Structure.objects.filter(
    entry__meta_data__value='icsd', 
    label='input')

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

# Convert oqmd structure to pymatgen structure
def pystruct(s):
    lattice = s.cell
    species = [i.species for i in s.atoms]
    coordinates = [i.coord for i in s.atoms]
    return python.Structure(lattice, species, coordinates)

print('update')