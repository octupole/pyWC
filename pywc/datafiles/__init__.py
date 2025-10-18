# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
"""
Datafiles module for pyWC test data.

Provides paths to test trajectories and structures for examples and testing.
"""

import os
import re

# Get the data directory path
_data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data')

# Test files for Willard-Chandler surface analysis
NPT_RUN_TPR = os.path.join(_data_dir, 'npt_run.tpr')
TRAJ_TEST_XTC = os.path.join(_data_dir, 'trj_test.xtc')
SELECTION_TXT = os.path.join(_data_dir, 'selection.txt')

# For backward compatibility with old examples (if needed)
# These will be the same files but with legacy names
MICELLE_PDB = NPT_RUN_TPR  # Placeholder - will need actual micelle file
WATER_GRO = NPT_RUN_TPR    # Placeholder - will need actual water file
WATER_XTC = TRAJ_TEST_XTC  # Using test trajectory


class Data(object):
    """A class for storing/accessing configurations, trajectories, topologies"""

    @staticmethod
    def sigeps(data, input_type):
        nm2angs = 10.0
        a, b = float(data[5]), float(data[6])
        sigma = 0
        if input_type == 'c6c12':
            c6, c12 = a, b
            if (c6 > 0.0):
                sigma = (c12 / c6)**(1. / 6.)
        else:
            sigma = a

        return sigma * nm2angs

    def __init__(self):
        self._label = list()
        self.label = list()
        self.file = dict()
        self.type = dict()
        self.format = dict()
        self.description = dict()

    def add(self, label, filetype, fileformat, desc):
        self._label.append(label)
        if label[0] != '_':
            self.label.append(label)
        self.file[label] = globals()[label]
        file = self.file[label]
        self.type[file] = filetype
        self.type[label] = filetype
        self.format[file] = fileformat
        self.format[label] = fileformat
        self.description[file] = desc
        self.description[label] = desc

    def vdwradii(self, filename):
        """Get van der Waals radii from a topology file"""
        if self.type.get(filename) == 'topol' and self.format.get(filename) == 'GMX':
            return self._vdwradii_gmx(filename)
        # Return empty dict if not a GMX topology file
        return {}

    def _vdwradii_gmx(self, filename):
        """Extract VDW radii from GROMACS topology file"""
        with open(filename) as f:
            input_type = 'sigeps'
            content = f.read()
            if re.match('.*name.*c6 *c12.*', content.replace('\n', ' ')):
                input_type = 'c6c12'
            f.seek(0)
            scan = False
            radii = dict()
            for line in f:
                if (scan and re.match(r'^ *\[', line)):
                    return radii
                if (scan):
                    try:
                        data = (line.split(";")[0]).split()
                        atom = data[0]
                        radii[atom] = 0.5 * self.sigeps(data, input_type)
                    except IndexError:
                        pass
                if (re.match(r'^ *\[ *atomtypes *\]', line)):
                    scan = True
        return radii

    def _generate_data_property(self, name):
        labels = []
        for label in self.type.keys():
            if self.type[label] == name:
                labels.append(label)
        return list(set(labels) & set(self.label))

    @property
    def config(self):
        return self._generate_data_property('config')

    @property
    def topol(self):
        return self._generate_data_property('topol')

    @property
    def traj(self):
        return self._generate_data_property('traj')


# Instantiate the pywc_data object
pywc_data = Data()

__all__ = [
    'NPT_RUN_TPR',
    'TRAJ_TEST_XTC',
    'SELECTION_TXT',
    'MICELLE_PDB',
    'WATER_GRO',
    'WATER_XTC',
    'pywc_data',
]
