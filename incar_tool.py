#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 13:03:53 2019

@author: mitchell koerner
@email:  mfkoerner@gmail.com
"""





class Incar():
    """
    Docstring
    """
    def __init__(self, system = None, read_fp = None):
        self.data = {}
        if not read_fp is None:
            self.read(read_fp)
        if system is None:
            self.system = 'unknown'
        else:
            self.system = system

    def read(self, filepath = 'INCAR', set_system = False):
        """
        reads INCAR to data
        """
        with open(filepath, 'r') as f:









