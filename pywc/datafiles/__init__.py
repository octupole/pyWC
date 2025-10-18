"""
Datafiles module for pyWC test data.

Provides paths to test trajectories and structures for examples and testing.
"""

import os

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

__all__ = [
    'NPT_RUN_TPR',
    'TRAJ_TEST_XTC',
    'SELECTION_TXT',
    'MICELLE_PDB',
    'WATER_GRO',
    'WATER_XTC',
]
