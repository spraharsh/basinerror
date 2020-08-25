#!/usr/bin/env python3
"""
Module for functions that determine whether the minimum is the same or not.
This has been adopted from CheckSameMinimum in basinvolume.

Assumes that the particle isn't a rattler
"""
from types import SimpleNamespace
import numpy as np




class CheckSameMinimum:
    """
    @brief      Container and methods for identifying minima

    @details    Container and methods for identifying minima
    """

    def __init__(self, ctol, dim, boxl=-1,
                 minimalist_max_len=1, minima_database_location=None, update_database=True):
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
        self.ctol = ctol
        self.dim = dim
        self.boxl = boxl

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

    def add_minimum(self, final_minimum, initial_coords,
                    failed_quench, part_ind=0):
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
        self.initial_coords_list.append(initial_coords)
        if failed_quench:
            # assign -1 to failed quench
            self.orderparamlist.append(-1)
        else:
            minima_asssigned = False
            for i, minimum in enumerate(self.minimalist):
                aligned_final_minimum = self.align_structures(minimum,
                                                              final_minimum,
                                                              part_ind)

                if self.check_same_structure(aligned_final_minimum, minimum):
                    self.orderparamlist.append(i)
                    minima_asssigned = True
                    break
            # assign new minimum if possible
            if (minima_asssigned is False):
                if len(self.minimalist) < self.minimalist_l:
                    self.minimalist.append(final_minimum)
                    self.orderparamlist.append(len(self.minimalist)-1)
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
        return boxed_coords.reshape(len(coords) // self.dim,
                                    self.dim)

    def check_same_structure(self, aligned_minimum, correct_minimum):
        """ checks whether the structure corresponding to
            the coordinates of the minima is the same
        Parameters
        ----------
        aligned_minimum: array
            aligned coordinates of the minimum to be compared
        correct_minimum: array
            coordinates of the correct minimum
        """
        dist_vec_array = (aligned_minimum - correct_minimum)
        # d_array[i] represents the distance between
        # the positions of the particle in the aligned minimum
        # and the correct minimum
        d_array = np.array(list(map(np.linalg.norm, dist_vec_array)))
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

    def align_structures(self, correct_minimum, obtained_minimum,
                         part_ind):
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
        driftless_min = map(lambda x: x-drift, obtained_minimum)
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
        file_fun = (lambda prefix, name : prefix +
                    '/' + name + maxminstr + '.npy')
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
            
        

    @staticmethod
    def load_map(foldpath, max_minima_l):
        """ loads data from file location

        

        Parameters
        ----------
        foldpath: path
            folder path
        max_minima_l: int
            max length of minimalist
        """
        maxminstr = str(max_minima_l)
        file_fun = (lambda prefix, name : prefix +
                    '/' + name  + maxminstr + '.npy')
        res = SimpleNamespace()
        res.initial_coords = np.load(file_fun(foldpath, 'initial_coords'))
        res.order_params = np.load(file_fun(foldpath, '/order_params'))
        res.minimalist = np.load(file_fun(foldpath, '/minimalist'))
        return res 


if __name__ == "__main__":
    th = CheckSameMinimum(1, 1)
    a = np.array([1.1, 2, 3])
    b = np.array([1, 2, 3])
    th.minimalist = [a]
    
