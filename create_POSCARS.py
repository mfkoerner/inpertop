#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 14:31:35 2019

@author: mitchell koerner
@email:  mfkoerner@gmail.com
"""

import get_structures as gs
from filter_structures import final_list
import os
from qmpy import io

print(final_list)
bonuses = [gs.InversePerovskiteBonuses(i) for i in final_list]

# make base directories
# for bonus in bonuses:
#     os.mkdir('/home/oqmd/pod/mkoerner/inpertop_data/icsd/I{}'.format(bonus.idstr))

# create POSCAR files
for bonus in bonuses:
    io.write(bonus.structure,
        filename='/home/oqmd/pod/mkoerner/inpertop_data/icsd/I{}/POSCAR'.format(bonus.idstr),
        comments = bonus.label
        )