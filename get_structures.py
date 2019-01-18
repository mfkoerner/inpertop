#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 18:05:05 2019

@author: mitchell koerner
@email:  mfkoerner@gmail.com
"""


import qmpy         #quantum materials database interface
import numpy
import os
import sys
from django.db.models import F
import importlib as imp

from qmpy import Structure

def upgs():
    import gs
    imp.reload(gs)
    from gs import *
    return()

allstruct = Structure.objects.filter(entry__meta_data__value__contains='perovskite', label='input')
icsd = Structure.objects.filter(
    entry__meta_data__value='icsd', 
    label='input')
unicsd53 = Structure.objects.filter(
    entry__meta_data__value='icsd', 
    label='input',
    natoms = 5,
    ntypes = 3,
    entry__id = F('entry__duplicate_of__id')
    )
