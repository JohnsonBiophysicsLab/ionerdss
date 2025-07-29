"""
This file includes code adapted from the 'pointgroup' package,
originally authored by Abel Carreras (https://github.com/abelcarreras/pointgroup),
and is licensed under the MIT License:

The MIT License (MIT)

Copyright (c) 2023 Efrem Bernuz and Abel Carreras

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""
import numpy as np
from . import utils
from .grid import get_cubed_sphere_grid_points
from ionerdss.utils.coords import get_perpendicular_vector
from ionerdss.utils.rotations import Rotation, rotation_matrix
from ionerdss.utils import inertia_tensors, angles

class PointGroup:
    """
    Point group main class. Note that we assume that center of mass is
    already put at (0,0,0).
    """

    def __init__(self,
                 positions,  # binding site positions
                 symbols,  # binding site symbols
                 tolerance_eig=1e-2,  # inertia tensor precision
                 tolerance_ang=4  # angular tolerance in degrees
                 ):

        self._tolerance_eig = tolerance_eig
        self._tolerance_ang = np.deg2rad(tolerance_ang)
        self._symbols = symbols
        self._cent_coord = np.array(positions)

        self._ref_orientation = np.identity(3)

        # determine inertia tensor
        inertia_tensor = inertia_tensors.get_inertia_tensor(self._cent_coord)
        eigenvalues, eigenvectors = np.linalg.eigh(inertia_tensor)
        self._eigenvalues = eigenvalues
        self._eigenvectors = eigenvectors.T

        # initialize variables
        self._schoenflies_symbol = ''
        self._max_order = 1

        eig_degeneracy = inertia_tensors.get_degeneracy(self._eigenvalues, self._tolerance_eig)

        # Linear groups
        if np.min(abs(self._eigenvalues)) < self._tolerance_eig:
            self._lineal()

        # Asymmetric group
        elif eig_degeneracy == 1:
            self._asymmetric()

        # Symmetric group
        elif eig_degeneracy == 2:
            self._symmetric()

        # Spherical group
        elif eig_degeneracy == 3:
            self._spherical()

        else:
            raise ValueError('Group type error')

    def _rename_point_group(self, pointgroup):
        """
        renames isomorphic point groups
        Ref) http://dms.library.utm.my:8080/vital/access/manager/Repository/vital:110227;jsessionid=8D0CBD4C6AD10FF6BEA7F6839D6E44BB?site_name=GlobalView

        :param pointgroup
        :return: renamed pointgroup
        """
        # C1h = Cs = C1v = S1
        if pointgroup in ['C1h', 'Cs', 'C1v', 'S1']:
            return 'Cs'
        # S2 = Ci
        elif pointgroup in ['S2', 'Ci']:
            return 'Ci'
        # C2 = D1
        elif pointgroup in ['C2', 'D1']:
            return 'C2'
        # C2h = D1d
        elif pointgroup in ['C2h', 'D1d']:
            return 'C2h'
        # C2v = D1h
        elif pointgroup in ['C2v', 'D1h']:
            return 'C2v'
        # C3h = S3 
        elif pointgroup in ['C3h', 'S3']:
            return 'C3h'
        else:
            return pointgroup

    def get_point_group(self):
        """
        get the point symmetry group symbol

        :return: the point symmetry group symbol
        """
        return self._rename_point_group(self._schoenflies_symbol)

    def get_standard_coordinates(self):
        """
        get the coordinates centered in the center of mass and
        oriented along principal axis of inertia

        :return: the coordinates
        """
        return self._cent_coord.tolist()

    def get_principal_axis_of_inertia(self):
        """
        get the principal axis of inertia in rows in increasing order of momenta of inertia

        :return: the principal axis of inertia
        """
        return self._eigenvectors.tolist()

    def get_principal_moments_of_inertia(self):
        """
        get the principal momenta of inertia in increasing order

        :return: list of momenta of inertia
        """
        return self._eigenvalues.tolist()

    # internal methods
    def _lineal(self):

        # set orientation
        idx = np.argmin(self._eigenvalues)
        main_axis = self._eigenvectors[idx]
        p_axis = get_perpendicular_vector(main_axis)
        self._set_orientation(main_axis, p_axis)

        # not considering reflection / inversion
        self._schoenflies_symbol = 'Cinfv'

    def _asymmetric(self):

        self._set_orientation(self._eigenvectors[2], self._eigenvectors[1])

        n_axis_c2 = 0
        main_axis = [1, 0, 0]
        for axis in np.identity(3):
            c2 = Rotation(axis, order=2)
            if self._check_op(c2, tol_factor=0.0):
                n_axis_c2 += 1
                main_axis = axis

        self._max_order = 2

        if n_axis_c2 == 0:
            self._max_order = 0
            self._no_rot_axis()
        else:
            self._cyclic(main_axis)

    def _symmetric(self):
        """
        handle cyclic groups Cn

        :return:
        """
        idx = utils.get_non_degenerated(self._eigenvalues, self._tolerance_eig)
        main_axis = self._eigenvectors[idx]

        self._max_order = self._get_axis_rot_order(main_axis, n_max=9)

        self._cyclic(main_axis) 

    def _spherical(self):
        """
        Handle spherical groups (T, O, I)

        :return:
        """

        main_axis = None
        while main_axis is None:
            for axis in get_cubed_sphere_grid_points(self._tolerance_ang):
                c5 = Rotation(axis, order=5)
                c4 = Rotation(axis, order=4)
                c3 = Rotation(axis, order=3)

                if self._check_op(c5, tol_factor=utils.magic_formula(5)):
                    self._schoenflies_symbol = "I"
                    main_axis = axis
                    self._max_order = 5
                    break
                elif self._check_op(c4, tol_factor=utils.magic_formula(4)):
                    self._schoenflies_symbol = "O"
                    main_axis = axis
                    self._max_order = 4
                    break
                elif self._check_op(c3, tol_factor=utils.magic_formula(3)):
                    self._schoenflies_symbol = "T"
                    main_axis = axis
                    self._max_order = 3

            if main_axis is None:
                print('increase tolerance')
                self._tolerance_ang *= 1.01

        p_axis_base = get_perpendicular_vector(main_axis)

        # I
        if self._schoenflies_symbol == 'I':
            def determine_orientation_I(main_axis):
                r_matrix = rotation_matrix(p_axis_base, np.arcsin((np.sqrt(5)+1)/(2*np.sqrt(3))))
                axis = np.dot(main_axis, r_matrix.T)

                # set molecule orientation in I
                for angle in np.arange(0, 2*np.pi / self._max_order+self._tolerance_ang, self._tolerance_ang):
                    rot_matrix = rotation_matrix(main_axis, angle)

                    c5_axis = np.dot(axis, rot_matrix.T)
                    c5 = Rotation(c5_axis, order=5)

                    if self._check_op(c5, tol_factor=utils.magic_formula(5)*np.sqrt(2)):
                        t_axis = np.dot(main_axis, rotation_matrix(p_axis_base, np.pi/2).T)
                        return np.dot(t_axis, rot_matrix.T)

                raise ValueError('Error orientation I group')

            p_axis = determine_orientation_I(main_axis)
            self._set_orientation(main_axis, p_axis)

        # O
        if self._schoenflies_symbol == 'O':

            # set molecule orientation in O
            def determine_orientation_O(main_axis):
                r_matrix = rotation_matrix(p_axis_base, np.pi/2)
                axis = np.dot(main_axis, r_matrix.T)

                for angle in np.arange(0, 2*np.pi / self._max_order+self._tolerance_ang, self._tolerance_ang):
                    rot_matrix = rotation_matrix(main_axis, angle)

                    c4_axis = np.dot(axis, rot_matrix.T)
                    c4 = Rotation(c4_axis, order=4)

                    if self._check_op(c4, tol_factor=utils.magic_formula(4)*np.sqrt(2)):
                        return axis

                raise ValueError('Error orientation O group')

            p_axis = determine_orientation_O(main_axis)
            self._set_orientation(main_axis, p_axis)

        # T
        if self._schoenflies_symbol == 'T':

            # set molecule orientation in T
            def determine_orientation_T(main_axis):
                r_matrix = rotation_matrix(p_axis_base, -np.arccos(-1/3))
                axis = np.dot(main_axis, r_matrix.T)

                for angle in np.arange(0, 2*np.pi / self._max_order + self._tolerance_ang, self._tolerance_ang):
                    rot_matrix = rotation_matrix(main_axis, angle)

                    c3_axis = np.dot(axis, rot_matrix.T)
                    c3 = Rotation(c3_axis, order=3)

                    if self._check_op(c3, tol_factor=utils.magic_formula(3)*np.sqrt(2)):
                        t_axis = np.dot(main_axis, rotation_matrix(p_axis_base, np.pi/2).T)
                        return np.dot(t_axis, rot_matrix.T)

                raise ValueError('Error orientation T group')
            
            p_axis = determine_orientation_T(main_axis)
            self._set_orientation(main_axis, p_axis)

    def _no_rot_axis(self):
        self._schoenflies_symbol = 'C1'
        return

    def _cyclic(self, main_axis):
        self._schoenflies_symbol = "C{}".format(self._max_order)
        return

    def _get_axis_rot_order(self, axis, n_max):
        """
        Get rotation order for a given axis

        :param axis: the axis
        :param n_max: maximum order to scan
        :return: order
        """

        def max_rotation_order(tolerance):
            for i in range(2, 100):
                if 2*np.pi / (i * (i - 1)) <= tolerance:
                    return i-1

        n_max = np.min([max_rotation_order(self._tolerance_ang), n_max])

        for i in range(n_max, 1, -1):
            Cn = Rotation(axis, order=i)
            if self._check_op(Cn):
                return i
        return 1

    def _check_op(self, operation, print_data=False, tol_factor=1.0):
        """
        check if operation exists

        :param operation: operation orbject
        :return: True or False
        """
        sym_matrix = operation.get_matrix()
        error_abs_rad = angles.absolute_error_to_angle(self._tolerance_eig, points=self._cent_coord)

        op_coordinates = np.dot(self._cent_coord, sym_matrix.T)
        for idx, op_coord in enumerate(op_coordinates):

            difference_rad = utils.normalized_radius_difference(op_coord, self._cent_coord, self._tolerance_eig)
            difference_ang = angles.angles_between_vector_and_vectors(op_coord, self._cent_coord, self._tolerance_eig)

            def check_diff(diff, diff2):
                for idx_2, (d1, d2) in enumerate(zip(diff, diff2)):
                    if self._symbols[idx_2] != self._symbols[idx]:
                        continue
                    # d_r = np.linalg.norm([d1, d2])
                    tolerance_total = self._tolerance_ang * tol_factor + error_abs_rad[idx_2]
                    if d1 < tolerance_total and d2 < tolerance_total:
                        return True
                return False

            if not check_diff(difference_ang, difference_rad):
                return False

            if print_data:
                print('Continue', idx)

        if print_data:
            print('Found!')
        return True

    def _set_orientation(self, main_axis, p_axis):
        """
        set molecular orientation along main_axis (x) and p_axis (y).

        :param main_axis: principal orientation axis (must be unitary)
        :param p_axis: secondary axis perpendicular to principal (must be unitary)
        :return:
        """

        assert np.linalg.norm(main_axis) > 1e-1
        assert np.linalg.norm(p_axis) > 1e-1

        orientation = np.array([main_axis, p_axis, np.cross(main_axis, p_axis)])
        self._cent_coord = np.dot(self._cent_coord, orientation.T)
        self._ref_orientation = np.dot(self._ref_orientation, orientation.T)
