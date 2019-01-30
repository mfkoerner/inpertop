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

bonuses = [gs.InversePerovskiteBonuses(i) for i in final_list]

for bonus in bonuses:
    os.mkdir('/home/oqmd/pod/mkoerner/inpertop_data/icsd/I{}'.format(bonus.idstr))
    