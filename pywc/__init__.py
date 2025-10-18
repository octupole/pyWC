# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
"""
Minimal pyWC namespace exposing the Willard–Chandler surface machinery.
"""

from .willard_chandler import WillardChandler
from . import utilities, datafiles
from .version import __version__

__all__ = ["WillardChandler", "utilities", "datafiles", "__version__"]

try:  # optional mdtraj tweak remains for compatibility
    from .patches import patchMDTRAJ_ReplacementTables
    patchMDTRAJ_ReplacementTables()
except Exception:
    pass
