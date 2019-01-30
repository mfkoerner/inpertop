#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 14:33:43 2019

@author: mitchell koerner
@email:  mfkoerner@gmail.com
"""

import get_structures as gs
from qmpy.materials.element import Element


# Useful objects
El_O = Element.get('O')
El_F = Element.get('F')

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
nonox = [ i for i in un53icsd221 if not 'O3' in i.name ]
# Knock out all with 3 Fluorine
nonoxF = [ i for i in nonox if not 'F3' in i.name ]

#Filter out atoms with most electronegativity in X site
invper_max = [i for i in nonoxF if gs.is_inv(i, 'max')]

# filter out F block electrons
invper_nof = [i for i in invper_max if len(gs.PARTIAL_F.intersection(set(i.elements))) == 0]

final_list = invper_nof