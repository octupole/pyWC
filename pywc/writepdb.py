# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
from __future__ import print_function
import numpy as np
import MDAnalysis


def _writepdb(interface,
              filename='layers.pdb',
              centered='no',
              group='all',
              multiframe=True,
              tempfactors=None):
    """

    Write the frame to a pdb file: please use :func:`Interface.writepdb`.

    Tests:

    >>> import pywc  # doctest: +SKIP
    >>> import MDAnalysis as mda  # doctest: +SKIP
    >>> from pywc.datafiles import NPT_RUN_TPR, TRAJ_TEST_XTC, SELECTION_TXT  # doctest: +SKIP
    >>> u = mda.Universe(NPT_RUN_TPR, TRAJ_TEST_XTC)  # doctest: +SKIP
    >>> with open(SELECTION_TXT) as f:  # doctest: +SKIP
    ...     selection = f.read().strip()  # doctest: +SKIP
    >>> g = u.select_atoms(selection)  # doctest: +SKIP
    >>> wc = pywc.WillardChandler(u, group=g, alpha=3.0, mesh=2.5)  # doctest: +SKIP
    >>> wc.writepdb('test.pdb', centered='no', multiframe=False)  # doctest: +SKIP
    >>> u2 = mda.Universe('test.pdb')  # doctest: +SKIP
    >>> print(f"Atoms written: {len(u2.atoms)}")  # doctest: +SKIP
    Atoms written: 105066
    >>>
    >>> # Skip remaining GITIM tests  # doctest: +SKIP
    >>> for centering in ['no', 'middle', 'origin']:  # doctest: +SKIP
    ... 	name='gitim.water.planar.'+centering+'.pdb'
    ... 	inter.writepdb(name,centered=centering,multiframe=False)
    ... 	u2 = mda.Universe(name)
    ... 	print(u2.atoms[0].position)
    [28.62  2.01 11.37]
    [28.62  2.01 62.88]
    [ 28.62   2.01 -12.12]
    >>>
    >>> print ('water itim')
    water itim
    >>> g = u.select_atoms('name OW')
    >>> inter = pywc.ITIM(u,group=g)  # doctest: +SKIP
    >>> for centering in ['no', 'middle', 'origin']:
    ... 	name='itim.water.'+centering+'.pdb'
    ... 	inter.writepdb(name,centered=centering,multiframe=False)
    ... 	u2 = mda.Universe(name)
    ... 	print(u2.atoms[0].position)
    [28.62  2.01 11.37]
    [28.62  2.01 62.88]
    [ 28.62   2.01 -12.12]
    """

    if isinstance(group, interface.universe.atoms.__class__):
        interface.group = group
    else:
        interface.group = interface.universe.atoms

    # we back up the original tempfactors and assign the new ones
    if tempfactors is not None:
        _tf = interface.group.atoms.tempfactors.copy()
        interface.group.tempfactors = tempfactors[:]

    temp_pos = np.copy(interface.universe.atoms.positions)
    move_to_origin = False
    options = {
        'no': False,
        False: False,
        'middle': True,
        'origin': True,
        True: True
    }
    if centered == 'origin':
        move_to_origin = True

    if options[centered] != interface.do_center:
        # i.e. we don't have already what we want ...
        if interface.do_center == False:  # we need to center
            interface.center(planar_to_origin=move_to_origin)
        else:  # we need to put back the original positions
            try:
                # original_positions are (must) always be defined
                interface.universe.atoms.positions = interface.original_positions
            except:
                raise AttributeError
    try:
        # it exists already, let's add information about the box, as
        # MDAnalysis forgets to do so for successive frames. A bugfix
        # should be on the way for the next version...
        interface.PDB[filename].CRYST1(
            interface.PDB[filename].convert_dimensions_to_unitcell(
                interface.universe.trajectory.ts))
    except:
        bondvalue = None
        interface.PDB[filename] = MDAnalysis.Writer(
            filename,
            multiframe=multiframe,
            n_atoms=interface.group.atoms.n_atoms,
            bonds=bondvalue)
    interface.PDB[filename].write(interface.group.atoms)
    interface.PDB[filename].pdbfile.flush()
    interface.universe.atoms.positions = np.copy(temp_pos)

    # we copy back the original tempfactors
    if tempfactors is not None:
        interface.group.tempfactors = _tf[:]
