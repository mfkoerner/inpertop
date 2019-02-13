#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 15:41:45 2019

@author: mitchell koerner
@email:  mfkoerner@gmail.com
"""

import numpy as np
import pandas as pd
import time
import os
from os.path import join as j

def datestring():
    now = time.localtime()
    out = '{:04d}_{:02d}_{:02d}'.format(now.tm_year, now.tm_mon, now.tm_mday)
    return out

def datetimestring():
    now = time.localtime()
    out = '{:04d}_{:02d}_{:02d}_{:02d}_{:02d}'.format(
        now.tm_year, now.tm_mon, now.tm_mday, now.tm_hour, now.tm_min
        )

class interface():
    """
    Documentation
    I'd higly recommend using an absolute filepath!
    """
    def __init__(self, filepath):
        self.filepath = filepath
        self.filename = os.path.basename(self.filepath)
        self.directory = os.path.dirname(self.filepath)
        self._read()
        self.recentDF = self.DF.copy()

    def _read(self):
        """
        loads csv from filepath to self.DF dataframe
        """
        self.DF = pd.read_csv(self.filepath, header = 0, index_col = 0)
        self.origDF = self.DF.copy()

    def _update_recent(self):
        self.recentDF = self.DF.copy()

    def write(self):
        """
        writes csv from dataframe to filepath
        """
        self.DF.to_csv(self.filepath)

    def backup(self, filepath = 'auto'):
        """
        backs up original dataframe to a backup file
        """
        if filepath == 'auto':
            backup_filepath = '{}/backups/{}-{}'.format(
                self.directory, datetimestring(), self.filename
                )
        else:
            backup_filepath = filepath
        self.origDF.to_csv(backup_filepath)  # note this differs from self.filepath

    def update_entry(self, entry, column, value):
        """
        updates a single entry (or multiple entries) in the dataframe
        """
        self.DF.loc[entry][column] = value

    def get_entry(self, entry, column = None):
        """
        gets a single entry (or multiple entries) and can be limited to certain columns
        """
        if column is None:
            return self.DF.loc[entry]
        else:
            return self.DF.loc[entry][column]

    def record_changes(self, directory = '.', update = True):
        """
        records changes to dataframe in a specified directory
        """
        updates = self.DF[self.DF != self.recentDF].dropna(how='all', inplace = False)
        filename = 'changes-{}-{}'.format(datetimestring(), self.filename)
        change_filepath = j(directory, filename)
        updates.to_csv(change_filepath)
        if update:
            self._update_recent()














