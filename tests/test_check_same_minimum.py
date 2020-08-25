"""
Tests for checksameminimum
"""
import pytest
import numpy as np
from ..checksameminimum import CheckSameMinimum
import os


def test_align():
    # print("hello world")
    specified_min = np.array([[1, 2, 3],
                              [3, 4, 5]])
    ctol = 0.1
    dim = 3
    min_check = CheckSameMinimum(ctol, dim)
    min_check.minimalist = [specified_min]
    dummy_quench_result = specified_min + np.array([[0.1, 0.2, 0.3],
                                                    [0, 0, 0]])
    aligned_quench = min_check.align_structures(specified_min, dummy_quench_result, 0)
    assert(np.all(aligned_quench == np.array([[1., 2., 3.],
                                              [2.9, 3.8, 4.7]])
                  )
           )


def test_check_same_structure_no_period():
    specified_min = np.array([[1, 2, 3],
                              [2, 3, 4]])
    ctol_a = 0.5 - 0.1
    ctol_b = 0.5 + 0.1
    dim = 3
    min_check_a = CheckSameMinimum(ctol_a, dim)
    min_check_b = CheckSameMinimum(ctol_b, dim)
    min_check_a.minimalist = [specified_min]
    min_check_b.minimalist = [specified_min]
    dummy_quench_result = specified_min + np.array([[0.0, 0.3, 0.4],
                                                    [0, 0, 0]])
    assert(min_check_a.check_same_structure(specified_min,
                                            dummy_quench_result)
           is False)
    assert(min_check_b.check_same_structure(specified_min,
                                            dummy_quench_result)
           is True)

def test_check_save_load():
    """
       Check whether data can be saved and loaded
    """
    specified_min = np.array([1, 2, 3, 4, 5, 6])
    ctol_a = 0.5
    dim = 3
    maxminl = 1
    min_check = CheckSameMinimum(ctol_a,
                                 dim, boxl=2.0,
                                 minimalist_max_len=maxminl)
    min_check.minimalist.append(min_check.box_reshape_coords(specified_min))
    min_check.add_minimum(specified_min, [0, 0], failed_quench=False)
    min_check.add_minimum(specified_min, [0, 0], failed_quench=False)
    print(min_check.initial_coords_list)
    print(min_check.orderparamlist)
    foldpath = os.getcwd() + '/testdata'
    os.makedirs(fossldpath, exist_ok=True)
    min_check.dump_map(foldpath)
    res = min_check.load_map(foldpath, max_minima_l=1)
    assert(np.all(res.initial_coords == np.array(min_check.initial_coords_list)))
    print(np.all(res.order_params == np.array(min_check.orderparamlist)))
    print(np.all(res.minimalist == np.array(min_check.minimalist)))



def test_add_minimum():
    specified_min = np.array([[1, 2, 3],
                              [2, 3, 4]])
    ctol = 0.4
    dim = 3
    min_check = CheckSameMinimum(ctol, dim, minimalist_max_len=2, boxl=2.0)
    specified_min_flat = specified_min.flatten()
    min_check.add_minimum(specified_min_flat, 1, False)

    dummy_quench_result_a = specified_min + np.array([[0.0, 0.1, 0.2],
                                                      [0, 0, 0]])
    dummy_quench_result_b = specified_min + np.array([[0.0, 0.5, 0.5],
                                                      [0, 0.0, 0.0]])
    dummy_quench_result_c = specified_min + np.array([[0.0, 3.0, 0.0],
                                                      [0, 0.0, 0.0]])
    # check for putting in the first minima
    min_check.add_minimum(dummy_quench_result_a.flatten(), 2, False)
    # check for identifying a minima within tol
    min_check.add_minimum(dummy_quench_result_b.flatten(), 3, False)
    # check for identifying false quench
    min_check.add_minimum(dummy_quench_result_b.flatten(), 4, True)
    # check for identifying a minima outside tol but with
    # periodic boundary conditions
    min_check.add_minimum(dummy_quench_result_c.flatten(), 5, False)
    print(min_check.initial_coords_list)
    print(min_check.orderparamlist)
    print(min_check.minimalist)
    # check the minima are as expected onece put back in the box. periodic
    # boundary conditions
    print(min_check.minimalist[0])
    assert(np.all(min_check.minimalist[0] == np.array([[1., 0., 1.],
                                                       [0., 1., 0.]])))
    assert(np.all(min_check.minimalist[1] == np.array([[1., 0.5, 1.5],
                                                       [0., 1., 0.]])))
    # check whether the right number of minima have been identified
    assert(len(min_check.minimalist) == 2)

    # check whether the quenches have been identified correct
    assert(min_check.orderparamlist == [0, 0, 1, -1, -2])
