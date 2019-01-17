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

from qmpy import Structure

allstruct = Structure.objects.filter(entry__meta_data__value__contains='perovskite', label='input')