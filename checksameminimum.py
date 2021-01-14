#!/usr/bin/env python3
"""
Module for functions that determine whether the minimum is the same or not.
This has been adopted from CheckSameMinimum in basinvolume.

Assumes that the particle isn't a rattler
"""
from types import SimpleNamespace
import numpy as np

from basinvolume.utils import in_hull, origin_in_hull_2d
import logging

class CheckSameMinimum:
    """
    @brief      Container and methods for identifying minima

    @details    Container and methods for identifying minima
    """
    def __init__(self,
                 ctol,
                 dim,
                 boxl=-1,
                 minimalist_max_len=1,
                 minima_database_location=None,
                 update_database=True,
                 rattler_check=False,
                 potential=None,
                 hs_radii=None,
                 max_nrattlers=None):
        """
        
        
        Parameters
        ----------
        ctol : tolerance for the quenched minimum to be part of the same basin,
               should ideally be lower than the quench tolerance of the minimum
        dim  : dimension of the space. helps reshaping arrays appropriately

        boxl : length of box, negative value assumes infinite box, positive
               value implies periodic boundary conditions with those conditions
        minima_database_location: filepath
               if the location is specified then we want to use the minima in the database
        update_database: bool
               if this is true then we keep updating the database every run.
        Returns
        -------
        out : None
        """
        #  !hasattr convoluted way of figuring out what hs_radii does
        if ((potential == None or (not hasattr(hs_radii, 'shape'))) and rattler_check == True):
            raise Exception('can\'t do rattler check without potential details')
        self.potential = potential
        self.hs_radii = hs_radii
        self.rattler_check = rattler_check
        self.ctol = ctol
        self.dim = dim
        self.boxl = boxl
        self.max_nrattlers = max_nrattlers
        self.nfluidstates=0
        self.nrattlermin = 0
        # list of minima as identified by the coords. if the minima database is not defined
        # it is assumed to be empty. if minima database is specified and update database
        # is true, we will update the database when we dump the map

        try:
            self.minimalist = np.load(minima_database_location)
        except:
            self.minimalist = []

        #  We also want to store the initial coords in the container
        #  because we want to be able to order accordingly

        self.initial_coords_list = []
        self.minima_database_location = minima_database_location
        self.update_database = update_database

        # orderparamlist[i] = arg where i is the corresponding argument
        # in the initial coordinates where arg is the is the argument
        # of the corresponding minimum in the minimalist
        # this is important because it defines
        # our map between the initial coords to the
        # minima

        self.orderparamlist = []
        self.minimalist_l = minimalist_max_len

    def add_minimum(self,
                    final_minimum,
                    initial_coords,
                    failed_quench,
                    part_ind=0):
        """ function that checks if the minimum corresponds to any previous
            minima in the list and adds it to the new minima
        Parameters
        ----------
        final_minimum : final minimum after the quench
        initial_coords : initial coords which are mapped to this minimum
                         only saved to define the map between a list of these
                         and the orderparamlist. can just be representative
        failed_quench : if the quench has failed
                        save the order parameter as (-1)

        Returns
        -------
        out : None
        """

        # we want reasonable dimensions on the coordinates first
        # because we care about the coordinates of each particle for these
        # algorithms. we also impose periodic boundary conditions

        final_minimum = self.box_reshape_coords(final_minimum)
        if self.rattler_check:
            # returns rattlers in case there are any
            # and whether state is jammed (i.e not_fluid)
            not_fluid, rattlers = self._find_rattlers(final_minimum.flatten())
            if np.any(rattlers==0):
                self.nrattlermin +=1
        else:
            rattlers = None
            not_fluid = False
        self.initial_coords_list.append(initial_coords)
        if failed_quench:
            # assign -1 to failed quench
            self.orderparamlist.append(-1)
        else:
            minima_asssigned = False
            for i, minimum in enumerate(self.minimalist):
                aligned_final_minimum = self.align_structures(
                    minimum, final_minimum, part_ind)
                if self.check_same_structure(aligned_final_minimum, minimum, rattlers):
                    self.orderparamlist.append(i)
                    minima_asssigned = True
                    break
            # assign new minimum if possible
            if (minima_asssigned is False):
                if len(self.minimalist) < self.minimalist_l:
                    if not_fluid:
                        self.minimalist.append(final_minimum)
                    else:
                        self.fluid += 1
                        self.minimalist.append(-100)
                    self.orderparamlist.append(len(self.minimalist) - 1)
                else:
                    self.orderparamlist.append(-2)
        return None

    def box_reshape_coords(self, coords):
        """ Reshapes and makes sure that the coords are
             put in a box. This is useful to reorganize the flat
             pele arrays by particle and make sure we're not breaking
             periodicity

        Parameters
        ----------
        coords: array
            input coordinates
        Returns
        ---------
        coords but put in an array with periodic boundary conditions imposed
        """
        if (self.boxl >= 0):
            boxed_coords = coords % self.boxl
        else:
            boxed_coords = coords.copy()
        return boxed_coords.reshape(len(coords) // self.dim, self.dim)

    def check_same_structure(self, aligned_minimum, correct_minimum, rattlers=None):
        """ checks whether the structure corresponding to
            the coordinates of the minima is the same
        Parameters
        ----------
        aligned_minimum: array
            aligned coordinates of the minimum to be compared
        correct_minimum: array
            coordinates of the correct minimum
        """
        dist_vec_array = (aligned_minimum - correct_minimum) % self.boxl
        dist_vec_array = np.where(dist_vec_array>self.boxl/2, self.boxl-dist_vec_array, dist_vec_array)
        # if the displacements are too large flip the box
        # (basically we want to get the absolute minimum possible)
        # if np.any(dist_vec_array>self.boxl/2):
        #     dist_vec_array = self.boxl-dist_vec_array
        # d_array[i] represents the distance between
        # the positions of the particle in the aligned minimum
        # and the correct minimum
        d_array = np.array(list(map(np.linalg.norm, dist_vec_array)))*rattlers
        if (np.max(d_array) < self.ctol):
            return True
        return False

    def put_in_box(self, minimum):
        """Parameters
           ----------
           minimum: array
               array to put in box
        """
        if (self.boxl <= 0):
            return minimum
        else:
            return minimum % self.boxl

    def align_structures(self, correct_minimum, obtained_minimum, part_ind):
        """ aligns the coordinates of one particle from the obtained_minimum
            with the coordinates of the
            other minimum
        Parameters
        ----------
        correct_minimum: array
            minimum to be checked against
        obtained_minimum: array
            minimum we want to check, it is assumed that periodic boundary
            conditions are imposed
        part_ind: integer
            index of particle we want to check against
        Returns
        -------
        aligned coordinates of the new minimum
        """
        drift = obtained_minimum[part_ind] - correct_minimum[part_ind]
        driftless_min = map(lambda x: x - drift, obtained_minimum)
        return np.array(list(driftless_min))

    def dump_map(self, foldpath, quench_name=None):
        """saves the parameters in files in the corresponding folderpath

           file names are tied to the maximum number of minima possible
           since that should basically keep us together
           if update database is True then we want to update the minimalist
           with the new minima. 
        
        Parameters
        ----------
        self: type
            description
        foldpath: path
            folder containing information
        """

        maxminstr = str(self.minimalist_l)
        file_fun = (
            lambda prefix, name: prefix + '/' + name + maxminstr + '.npy')
        np.save(file_fun(foldpath, 'initial_coords'),
                np.array(self.initial_coords_list))
        np.save(file_fun(foldpath, '/order_params'),
                np.array(self.orderparamlist))
        if self.update_database is False:
            np.save(file_fun(foldpath, '/minimalist'),
                    np.array([a.flatten() for a in self.minimalist]))
        else:
            np.save(self.minima_database_location,
                    np.array([a.flatten() for a in self.minimalist]))


    def _find_rattlers(self, coords):
        """
        finish this, I need to remove the rattler and break.
        Also need to get compare to existing jammed_packing option
        :return:
        (Jammed state or not, Rattlers)
        copied johannes code from basinvolume

        """
        if self.dim < 4:
            zmin = self.dim + 1
        else:
            raise NotImplementedError
        if self.max_nrattlers == None:
            self.max_nrattlers = len(self.hs_radii)
        self.nparticles = len(self.hs_radii)/self.dim
        self.rattlers = np.zeros(len(self.hs_radii))
        check_again = range(len(self.hs_radii))
        neighbor_indicess, neighbor_distancess \
            = self.potential.getNeighbors(coords)
        nrattlers = 0
        if self.dim != 2:
            origin = np.zeros(self.dim)
        while len(check_again) > 0:
            check_inds = check_again
            check_again = set()
            if nrattlers > self.max_nrattlers:
                logging.warning("Too many rattlers. Discarding packing.")
                return False
            for atomi in check_inds:
                found_rattler = False
                no_neighbors = len(neighbor_indicess[atomi])
                if no_neighbors < zmin:
                    found_rattler = True
                    logging.debug("Particle {} is not isostatic."
                                            .format(atomi))
                else:
                    if self.dim == 2:
                        found_rattler = not origin_in_hull_2d(
                            neighbor_distancess[atomi])
                    else:
                        points = (np.asarray(neighbor_distancess[atomi])
                                  .reshape((-1, self.dim)))
                        found_rattler = not in_hull(origin, points)
                    if found_rattler:
                        logging.debug("Particle {} is not in "
                                                "contacts' convex hull."
                                                .format(atomi))
                self.rattlers[atomi] = False if found_rattler else True
                if found_rattler:
                    if atomi in check_again:
                        check_again.remove(atomi)
                    for atomj in neighbor_indicess[atomi]:
                        check_again.add(atomj)
                        i_in_j = neighbor_indicess[atomj].index(atomi)
                        del neighbor_indicess[atomj][i_in_j]
                        del neighbor_distancess[atomj][i_in_j]
                    nrattlers += 1
        # test that number of contacts is sufficient for bulk modulus to be positive,
        # see eq 4 in http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.109.095704
        # see eq 19 in arXiv:1406.1529
        print(nrattlers, 'nrattlers')
        no_stable = self.nparticles - nrattlers
        total_contacts = sum([len(neighbor_indices)
                              for neighbor_indices in neighbor_indicess])
        N_min = int(2 * (self.dim * (no_stable - 1) + 1))
        logging.debug("N_min: {} total_contacts: {}"
                      .format(N_min, total_contacts))
        logging.debug("Number of rattlers: {}".format(nrattlers))
        if nrattlers > self.max_nrattlers:
            logging.warning("Too many rattlers. Discarding packing.")
            return (False, self.rattlers)
        if total_contacts >= N_min:
            return (True, self.rattlers)
        else:
            logging.warning("Packing is not globally stable, N_min: {} "
                            "total_contacts: {}".format(N_min, total_contacts))
            return (False, self.rattlers)


    @staticmethod
    def load_map(foldpath, max_minima_l, minima_database_path=None):
        """ loads data from file location
        Parameters
        ----------
        foldpath: path
            folder path
        max_minima_l: int
            max length of minimalist
        """
        maxminstr = str(max_minima_l)
        file_fun = (
            lambda prefix, name: prefix + '/' + name + maxminstr + '.npy')
        res = SimpleNamespace()
        res.initial_coords = np.load(file_fun(foldpath, 'initial_coords'))
        res.order_params = np.load(file_fun(foldpath, '/order_params'))
        if minima_database_path == None:
            res.minimalist = np.load(file_fun(foldpath, '/minimalist'))
        else:
            res.minimalist = np.load(minima_database_path)
        return res


if __name__ == "__main__":
    th = CheckSameMinimum(1, 1)
    a = np.array([1.1, 2, 3])
    b = np.array([1, 2, 3])
    th.minimalist = [a]
