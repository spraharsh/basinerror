"""
Contains helper functions for figuring out accuracy of our basins
"""

import numpy as np
from checksameminimum import CheckSameMinimum
from map_basin_steepest import QUENCH_FOLDER_NAME, MINIMA_DATABASE_NAME
from params import BASE_DIRECTORY
import os


CORRECT_MINIMA_FOLDER = 'correct_minima'



def percentage_same(minima_path_a, minima_path_b, foldnameInversePower, ctol=1e-3, minima_l=2000):
    """
    Checks the percentage of minima that are same in path a and path b, and returns the fraction of minima
    that are not the same.
    """
    abs_path_a = BASE_DIRECTORY + '/' + foldnameInversePower + '/' + minima_path_a
    abs_path_b = BASE_DIRECTORY + '/' + foldnameInversePower + '/' + minima_path_b
    minima_database_path = BASE_DIRECTORY + '/' + foldnameInversePower + '/' + MINIMA_DATABASE_NAME
    data_a = CheckSameMinimum.load_map(abs_path_a, max_minima_l=minima_l, minima_database_path=minima_database_path)
    data_b = CheckSameMinimum.load_map(abs_path_b, max_minima_l=minima_l, minima_database_path=minima_database_path)
    op_a_arr = data_a.order_params
    op_b_arr = data_b.order_params
    return (sum(op_a_arr!=op_b_arr))/len(op_b_arr) 


def correct_minima_check(minima_path_a, foldnameInversePower):
    """
    Checks the percentage of minima that are same in path a and path b
    """
    abs_path_correct = BASE_DIRECTORY + '/' + foldnameInversePower + '/' + CORRECT_MINIMA_FOLDER
    minima_database_path = BASE_DIRECTORY + '/' + foldnameInversePower + '/' + MINIMA_DATABASE_NAME
    if not os.path.isdir(abs_path_correct):
        exception_string = 'correct minima folder '+ abs_path_correct + ' does not exist. needs to be generated or set correctly'
        raise Exception(exception_string)
    return percentage_same(minima_path_a, CORRECT_MINIMA_FOLDER, foldnameInversePower)

    
if __name__=="__main__":
    foldnameInversePower = "ndim=2phi=0.9seed=0n_part=16r1=1.0r2=1.4rstd1=0.05rstd2=0.06999999999999999use_cell_lists=0power=2.5eps=1.0"
    old_correct_minima = 'scratch'
    percentage_same(CORRECT_MINIMA_FOLDER, old_correct_minima, foldnameInversePower)
    print(correct_minima_check(old_correct_minima, foldnameInversePower), 'correct minima')
