# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
from abc import ABCMeta, abstractmethod
import numpy as np
from .properties import _create_property
from .writepdb import _writepdb
from . import messages
from . import utilities
from scipy.spatial import cKDTree

# Import optimized centering implementation
try:
    from ._center_impl import center_wrapper, HAS_CENTER_FAST
except ImportError:
    # Fallback if the module is not available
    HAS_CENTER_FAST = False
    center_wrapper = None


class Interface(object):
    """ The Interface metaclass. Classes for interfacial determination
        (ITIM, GITIM,...) are derived from this one
    """
    __metaclass__ = ABCMeta

    directions_dict = {
        0: 'x',
        1: 'y',
        2: 'z',
        'x': 'x',
        'y': 'y',
        'z': 'z',
        'X': 'x',
        'Y': 'y',
        'Z:': 'z'
    }
    symmetry_dict = {
        'generic': 'generic',
        'cylindrical': 'cylindrical',
        'spherical': 'spherical',
        'planar': 'planar'
    }

    # main properties shared by all implementations of the class
    # When required=True is passed, the implementation of the class *must*
    # override the method when instantiating the class (i.e., before __init__)
    # By default required=False, and the name is set to None

    # interface *must* be created first.
    alpha, _alpha =\
        _create_property('alpha', "(float) real space cutoff")
    layers, _layers =\
        _create_property('layers', "AtomGroups of atoms in layers")
    analysis_group, _analysis_group =\
        _create_property('analysis_group', "(AtomGroup) the group, "
                         "the surface of which should be computed")
    cluster_cut, _cluster_cut =\
        _create_property('cluster_cut', "(real) cutoff for phase "
                         "identification")
    molecular, _molecular =\
        _create_property('molecular', "(bool) whether to compute "
                         "surface atoms or surface molecules")
    surfaces, _surfaces =\
        _create_property('surfaces', "Surfaces associated to the interface",
                         readonly=True)
    info, _info =\
        _create_property('info', "(bool) print additional information")

    multiproc, _multiproc =\
        _create_property('multiproc', "(bool) use parallel implementation")

    extra_cluster_groups, _extra_cluster_groups =\
        _create_property('extra_cluster_groups',
                         "(ndarray) additional cluster groups")
    radii_dict, _radii_dict =\
        _create_property('radii_dict', "(dict) custom atomic radii")

    max_layers, _max_layers =\
        _create_property('max_layers',
                         "(int) maximum number of layers to be identified")
    autoassign, _autoassign =\
        _create_property('autoassign',
                         "(bool) assign layers every time a frame changes")
    include_zero_radius, _include_zero_radius=\
        _create_property('include_zero_radius',
                         "(bool) include atoms with zero radius in the analysis (excluded by default)")
    cluster_threshold_density, _cluster_threshold_density =\
        _create_property('cluster_threshold_density',
                         "(float) threshold for the density-based filtering")
    surface_cluster_cut, _surface_cluster_cut =\
        _create_property('surface_cluster_cut',
                         "(float) distance cut to calculate the surface clusters")

    # TODO: does this belong here ?
    _interpolator, __interpolator =\
        _create_property('_interpolator', "(dict) custom atomic radii")

    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def _assign_layers(self):
        pass

    @property
    def atoms(self):
        if len(self._layers) == 0: return self._layers # an empty atom group
        return self._layers[:].sum()

    @property
    def method(self):
        return self.__class__.__name__

    def label_planar_sides(self):
        """ Assign to all layers a label (the beta tempfactor)
            that can be used in pdb files. Additionally, set
            the new layers and sides.
        """
        for uplow in [0, 1]:
            for nlayer, layer in enumerate(self._layers[uplow]):
                if layer is None:
                    self._layers[uplow][nlayer] = self.universe.atoms[:0]
                else:
                    self.label_group(
                        layer, beta=nlayer + 1.0, layer=nlayer + 1, side=uplow)

    def label_group(self,
                    group,
                    beta=None,
                    layer=None,
                    cluster=None,
                    side=None):
        if group is None:
            raise RuntimeError(
                'one of the groups, possibly a layer one, is None.' +
                ' Something is wrong...')
        if len(group) == 0:
            return
        if self.molecular is True:
            _group = group.residues.atoms
        else:
            _group = group

        if beta is not None:
            _group.tempfactors = float(beta)
        if layer is not None:
            _group.layers = layer
        if side is not None:
            _group.sides = side
        if cluster is not None:
            _group.clusters = cluster

    def _assign_symmetry(self, symmetry):
        if self.analysis_group is None:
            raise TypeError(messages.UNDEFINED_ANALYSIS_GROUP)
        if symmetry == 'guess':
            raise ValueError("symmetry 'guess' To be implemented")
        else:
            if not (symmetry in self.symmetry_dict):
                raise ValueError(messages.WRONG_DIRECTION)
            self.symmetry = symmetry

    def _generate_surface_clusters(self, group, cut):
        # at the moment, selects only the biggest cluster
        labels, counts, neighs = utilities.do_cluster_analysis_dbscan(
                group, cut,
                molecular=False)
        return group[np.where(labels == np.argmax(counts))[0]]

    def _define_cluster_group(self):
        self.universe.atoms.pack_into_box()
        self.cluster_group = self.universe.atoms[:0]  # empty
        if (self.cluster_cut is not None):
            cluster_cut = float(self.cluster_cut[0])
            # we start by adding the atoms in the smaller clusters
            # of the opposit phase, if extra_cluster_groups are provided
            if (self.extra_cluster_groups is not None):
                for extra in self.extra_cluster_groups:
                    x_labels, x_counts, _ = utilities.do_cluster_analysis_dbscan(
                        extra, cluster_cut, self.cluster_threshold_density,
                        self.molecular)
                    x_labels = np.array(x_labels)
                    x_label_max = np.argmax(x_counts)
                    x_ids_other = np.where(x_labels != x_label_max)[0]

                    self.cluster_group += extra[x_ids_other]

            # next, we add the atoms belonging to the main phase
            self.cluster_group += self.analysis_group

            # groups have been checked already in _sanity_checks()
            # self.cluster_group at this stage is composed of analysis_group +
            # the smaller clusters of the other phase
            labels, counts, neighbors = utilities.do_cluster_analysis_dbscan(
                self.cluster_group, cluster_cut,
                self.cluster_threshold_density, self.molecular)
            labels = np.array(labels)

            # counts is not necessarily ordered by size of cluster.
            sorting = np.argsort(counts,kind='stable')[::-1]
            # labels for atoms in each cluster starting from the largest
            unique_labels = np.sort(np.unique(labels[labels > -1]))
            # by default, all elements of the cluster_group are in
            # single-molecule/atom clusters. We will update them right after.
            self.label_group(self.cluster_group, cluster=-1)
            # we go in reverse order to let smaller labels (bigger clusters)
            # overwrite larger labels (smaller cluster) when the molecular
            # option is used.
            for el in unique_labels[::-1]:
                # select a label
                cond = np.where(labels == el)
                if self.molecular is True:
                    g_ = self.cluster_group[cond].residues.atoms
                else:
                    g_ = self.cluster_group[cond]
                # probably we need an example here, say:
                # counts = [ 61, 1230, 34, 0, ...  0 ,0 ]
                # labels = [ 0, 1, 2, 1, -1  ....  -1 ]
                # we have three clusters, of 61, 1230 and 34 atoms.
                # There are 61 labels '0'
                #         1230 labels '1'
                #           34 labels '2'
                #         the remaining are '-1'
                #
                # sorting = [1,0,2,3,....] i.e. the largest element is in
                #     (1230) position 1, the next (61) is in position 0, ...
                # Say, g_ is now the group with label '1' (the biggest cluster)
                # Using argwhere(sorting==1) returns exactly 0 -> the right
                # ordered label for the largest cluster.
                self.label_group(g_, cluster=np.argwhere(sorting == el)[0, 0])
            # now that labels are assigned for each of the clusters,
            # we can restric the cluster group to the largest cluster.


            if self.min_cluster_size is not None and self.n_clusters is not None:
                raise ValueError("min_cluster_size and n_clusters cannot both take a value different from None")
            if self.n_clusters is not None:
                self.cluster_group = self.universe.atoms[:0] # reset to the empty group
                self._unique_labels = unique_labels[:]
                for el in unique_labels[:self.n_clusters]:
                    self.cluster_group += self.universe.atoms[self.universe.atoms.clusters == el]
            elif self.min_cluster_size is not None:
                self.cluster_group = self.universe.atoms[:0] # reset to the empty group
                for el in unique_labels:
                    g_ = self.universe.atoms[self.universe.atoms.clusters==el]
                    if len(g_) >= self.min_cluster_size:
                        self.cluster_group += g_
                    else:
                        break # unique_labels are sorted by size

            else: # we still filter out molecules which do not belong to any cluster
                ids = np.where(labels != -1)[0]
                self.cluster_group = self.cluster_group[ids]

            self.n_neighbors = neighbors
        else:
            self.cluster_group = self.analysis_group
            self.label_group(self.cluster_group, cluster=0)

    def is_buried(self, pos):
        """ Checks wether an array of positions are located below
            the first interfacial layer """
        inter = self
        box = inter.universe.dimensions[:3]
        nonsurface = inter.cluster_group - inter.atoms[inter.atoms.layers == 1]
        # there are no inner atoms, distance is always > 0
        if len(nonsurface) == 0:
            return np.asarray([True] * len(pos))
        tree = cKDTree(nonsurface.positions, boxsize=box)
        neighs = tree.query_ball_point(pos, inter.alpha)
        condition = np.array([len(el) != 0 for el in neighs])
        return condition

    def reset_labels(self):
        """ Reset labels before interfacial analysis"""
        self.label_group(
            self.universe.atoms, beta=0.0, layer=-1, cluster=-1, side=-1)

    @staticmethod
    def _attempt_shift(group, _pos_group, direction, halfbox_shift, _dir):
        if _dir == 'x':
            utilities.centerbox(
                group.universe,
                x=_pos_group,
                center_direction=direction,
                halfbox_shift=halfbox_shift)
        if _dir == 'y':
            utilities.centerbox(
                group.universe,
                y=_pos_group,
                center_direction=direction,
                halfbox_shift=halfbox_shift)
        if _dir == 'z':
            utilities.centerbox(
                group.universe,
                z=_pos_group,
                center_direction=direction,
                halfbox_shift=halfbox_shift)

    def prepare_box(self):
        """ Before the analysis, pack every molecule into the box.
            Keep the original positions for latter use.
        """
        self.original_positions = np.copy(self.universe.atoms.positions[:])
        self.universe.atoms.pack_into_box()

    @staticmethod
    def _center(group, direction, halfbox_shift=False, centering_backend='cpu'):
        """
        Centers the liquid slab in the simulation box.

        The algorithm tries to avoid problems with the definition
        of the center of mass. First, a rough density profile
        (10 bins) is computed. Then, the group is shifted
        and reboxed until the bins at the box boundaries have a
        density lower than a threshold delta


        In ITIM, the system along the normal direction is always
        centered at 0 (halfbox_shift==True). To center to the middle
        of the box along all directions, set halfbox_shift=False

        """
        # Use optimized C++ implementation if available
        if HAS_CENTER_FAST and center_wrapper is not None:
            use_gpu = (centering_backend.lower() == 'gpu' or centering_backend.lower() == 'cupy')
            return center_wrapper(group, direction, halfbox_shift, force_python=False, use_gpu=use_gpu)

        # Fallback to pure Python implementation
        dim = group.universe.coord.dimensions
        total_shift = 0

        if not (direction in Interface.directions_dict):
            raise ValueError(messages.WRONG_DIRECTION)
        _dir = Interface.directions_dict[direction]
        _xyz = {
            'x': (0, utilities.get_x),
            'y': (1, utilities.get_y),
            'z': (2, utilities.get_z)
        }

        if _dir in _xyz.keys():
            direction = _xyz[_dir][0]
            _pos_group = (_xyz[_dir][1])(group)

        shift = dim[direction] / 100.

        _x = utilities.get_x(group.universe.atoms)
        _y = utilities.get_y(group.universe.atoms)
        _z = utilities.get_z(group.universe.atoms)

        _range = (0., dim[direction])
        if (halfbox_shift is True):
            _range = (-dim[direction] / 2., dim[direction] / 2.)

        histo, _ = np.histogram(
            _pos_group, bins=10, range=_range, density=True)

        max_val, min_val = np.amax(histo), np.amin(histo)
        # NOTE maybe allow user to set different values
        delta = min_val + (max_val - min_val) / 3.

        # let's first avoid crossing pbc with the liquid phase. This can fail:
        while (histo[0] > delta or histo[-1] > delta):
            total_shift += shift
            if total_shift >= dim[direction]:
                raise ValueError(messages.CENTERING_FAILURE)
            _pos_group += shift

            Interface._attempt_shift(group, _pos_group, direction,
                                     halfbox_shift, _dir)

            histo, _ = np.histogram(
                _pos_group, bins=10, range=_range, density=True)

        _center_ = np.average(_pos_group)
        if (halfbox_shift is False):
            box_half = dim[direction] / 2.
        else:
            box_half = 0.
        _pos = {'x': _x, 'y': _y, 'z': _z}
        if _dir in _pos.keys():
            _pos[_dir] += total_shift - _center_ + box_half
        # finally, we copy everything back
        group.universe.atoms.positions = np.column_stack((_x, _y, _z))

    @staticmethod
    def shift_positions_to_middle(universe, normal):
        box = universe.dimensions[normal]
        translation = [0, 0, 0]
        translation[normal] = box / 2.
        universe.atoms.positions += np.array(translation)
        universe.atoms.pack_into_box()

    def _shift_positions_to_middle(self):
        Interface.shift_positions_to_middle(self.universe, self.normal)

    @staticmethod
    def center_system(symmetry, group, direction, planar_to_origin=False, centering_backend='cpu'):
        if symmetry == 'planar':
            utilities.centerbox(group.universe, center_direction=direction)
            Interface._center(group, direction, halfbox_shift=True, centering_backend=centering_backend)
            utilities.centerbox(group.universe, center_direction=direction)
            if planar_to_origin is False:
                Interface.shift_positions_to_middle(group.universe, direction)
        else:
            for xyz in [0, 1, 2]:
                try:
                    Interface._center(group, xyz, halfbox_shift=False, centering_backend=centering_backend)
                except ValueError:
                    pass
            group.universe.atoms.pack_into_box()

    def center(self, planar_to_origin=False):
        centering_backend = getattr(self, 'centering_backend', 'cpu')
        Interface.center_system(
            self.symmetry,
            self.cluster_group,
            self.normal,
            planar_to_origin=planar_to_origin,
            centering_backend=centering_backend)

        self.centered_positions = np.copy(self.universe.atoms.positions[:])

    def writepdb(self,
                 filename='layers.pdb',
                 centered='no',
                 group='all',
                 multiframe=True,
                 tempfactors=None):
        """ Write the frame to a pdb file, marking the atoms belonging
            to the layers with different beta factors.

            :param str       filename:    the output file name
            :param str       centered:    'origin', 'middle', or 'no'
            :param AtomGroup group:       if 'all' is passed, use universe
            :param bool      multiframe:  append to pdb file if True
            :param ndarray   tempfactors: use this array as temp (beta) factors

            Example: save the positions (centering the interface in the cell)

            >>> import pywc  # doctest: +SKIP
            >>> import MDAnalysis as mda  # doctest: +SKIP
            >>> from pywc.datafiles import NPT_RUN_TPR, TRAJ_TEST_XTC, SELECTION_TXT  # doctest: +SKIP
            >>> u = mda.Universe(NPT_RUN_TPR, TRAJ_TEST_XTC)  # doctest: +SKIP
            >>> with open(SELECTION_TXT) as f:  # doctest: +SKIP
            ...     selection = f.read().strip()  # doctest: +SKIP
            >>> group = u.select_atoms(selection)  # doctest: +SKIP
            >>> wc = pywc.WillardChandler(u, group=group, alpha=3.0, mesh=2.5)  # doctest: +SKIP
            >>> wc.writepdb('surface.pdb', multiframe=False)  # doctest: +SKIP
        """

        _writepdb(
            self,
            filename=filename,
            centered=centered,
            group=group,
            multiframe=multiframe,
            tempfactors=tempfactors)

    @staticmethod
    def _():
        """
        This is a collection of basic tests to check
        that code is running -- no test on the correctness
        of the output is performed here.

        >>> # TEST:0 loading the module
        >>> import pywc  # doctest: +SKIP
        >>> import MDAnalysis as mda  # doctest: +SKIP

        >>> # TEST:1 basic WillardChandler functionality
        >>> from pywc.datafiles import NPT_RUN_TPR, TRAJ_TEST_XTC, SELECTION_TXT  # doctest: +SKIP
        >>> u = mda.Universe(NPT_RUN_TPR, TRAJ_TEST_XTC)  # doctest: +SKIP
        >>> with open(SELECTION_TXT) as f:  # doctest: +SKIP
        ...     selection = f.read().strip()  # doctest: +SKIP
        >>> group = u.select_atoms(selection)  # doctest: +SKIP
        >>> wc = pywc.WillardChandler(u, group=group, alpha=3.0, mesh=2.5)  # doctest: +SKIP
        >>> print(f"Surface area: {wc.surface_area:.2f}")  # doctest: +SKIP
        Surface area: 29527.69

        """
        pass

    @staticmethod
    def __():
        """
        Placeholder for advanced tests.

        >>> # Basic WillardChandler test
        >>> import pywc  # doctest: +SKIP
        >>> import MDAnalysis as mda  # doctest: +SKIP
        >>> from pywc.datafiles import NPT_RUN_TPR, TRAJ_TEST_XTC, SELECTION_TXT  # doctest: +SKIP
        >>> u = mda.Universe(NPT_RUN_TPR, TRAJ_TEST_XTC)  # doctest: +SKIP
        >>> with open(SELECTION_TXT) as f:  # doctest: +SKIP
        ...     selection = f.read().strip()  # doctest: +SKIP
        >>> group = u.select_atoms(selection)  # doctest: +SKIP
        >>> wc = pywc.WillardChandler(u, group=group, alpha=3.0, mesh=2.5)  # doctest: +SKIP

        """
        pass


#
