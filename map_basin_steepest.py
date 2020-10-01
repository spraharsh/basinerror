"""
Maps the basin using the steepest descent algorithm. The accuracy of the mapping depends on the step size.
Note that our goal is to return the minimum given by the trajectory of the solution
 to d x/ d t = - \grad{E} (check the truth of this statement)
"""

import numpy as np
import os
from basinerror import quench_cvode_opt, quench_steepest
from params import load_params, load_secondary_params
from pele.potentials import InversePower
from params import BASE_DIRECTORY
import matplotlib.pyplot as plt
from pele.optimize._quench import modifiedfire_cpp
from basinerror import quench_mixed_optimizer
from checksameminimum import CheckSameMinimum
QUENCH_FOLDER_NAME = 'mxopt'
MINIMA_DATABASE_NAME = 'minima_database.npy'


def map_binary_inversepower(quench,
                            foldername,
                            particle_coords,
                            coordarg,
                            tol=0.1):
    """
    Finds whether a point defined by particle_coord
    on the meshgrid correspond to a minimum or not for a 2d
    case.
    """
    foldpath = BASE_DIRECTORY + '/' + foldername
    # import params
    sysparams = load_params(foldpath)
    (hs_radii, initial_coords, box_length) = load_secondary_params(foldpath)
    assert (sysparams.ndim.value == 2)
    minimum_coords = np.loadtxt(foldpath + '/coords_of_minimum.txt',
                                delimiter=',')
    quench_coords = minimum_coords.copy()
    quench_coords[coordarg] = particle_coords[0]
    quench_coords[coordarg + 1] = particle_coords[1]

    # print(quench_coords, 'quench coords')
    # box length
    box_length = float(box_length)
    boxv = [box_length] * sysparams.ndim.value
    potential = InversePower(sysparams.power.value,
                             sysparams.eps.value,
                             use_cell_lists=False,
                             ndim=sysparams.ndim.value,
                             radii=hs_radii * 1.0,
                             boxvec=boxv)
    # ret = quench_mixed_optimizer(potential,
    #                              quench_coords,  # make sure right coords are being passed
    #                              # stepsize=0.001,
    #                              T=10,
    #                              # step=1,
    #                              nsteps=10000,
    #                              # conv_tol=1e-2,
    #                              tol=1e-4)
    # ret = quench_steepest(
    #     potential,
    #     quench_coords,  # make sure right coords are being passed
    #     nsteps=2000000,
    #     stepsize=5e-3,  # for steepest descent step size should be small
    #     tol=1e-4)
    ret = quench_cvode_opt(potential, quench_coords, tol=1e-4)
    # ret = modifiedfire_cpp(quench_coords, potential, tol=1e-4)
    print(ret.nfev)
    results = (ret.coords, ret.success, coordarg, ret.nfev, ret.nsteps, ret.nhev)
    return results


def construct_point_set_2d(foldername, nmesh, boxlscale, coordarg):
    """ Constructs a point set in 2d around a minimum for mapping the
        basin projection for a single particle.
     Parameters
    ----------
        foldername : name of the folder containing
                   1) the coordinates of the minimum
                      around which we constuct the grid
                   2) the boxlength in a parameter file

        nmesh : the number of mesh points in 2 directions
        boxlscale : scale of the box
    Returns
    -------
        out : a set of points over which we wish to map the basin
    """

    foldpath = BASE_DIRECTORY + '/' + foldername
    (hs_radii, initial_coords, box_length) = load_secondary_params(foldpath)

    # get necessary parameters from the folder
    box_length = float(box_length)
    print(box_length)
    minimum_coords = np.loadtxt(foldpath + '/coords_of_minimum.txt',
                                delimiter=',')

    # sets the length of the x and y range over which
    # the mesh is calculated
    xylength = boxlscale * box_length
    print(xylength)
    print(minimum_coords, 'minimum')
    print(minimum_coords[coordarg])

    #  This division in cases is because linspace doesn't
    #  give the right answer for n = 1
    if (nmesh != 1):
        # x_range = np.linspace(minimum_coords[coordarg] - xylength/2,
        #   minimum_coords[coordarg] + xylength/2, nmesh)
        # y_range = np.linspace(minimum_coords[coordarg+1] - xylength/2,
        #   minimum_coords[coordarg+1] + xylength/2, nmesh)

        # position minimum in the center
        centered_x_of_min = minimum_coords[coordarg] + box_length / 2
        centered_y_of_min = minimum_coords[coordarg + 1] + box_length / 2

        x_range = centered_x_of_min + np.linspace(0, box_length, nmesh)
        y_range = centered_y_of_min + np.linspace(0, box_length, nmesh)
    else:
        x_range = np.array([minimum_coords[coordarg]])
        y_range = np.array([minimum_coords[coordarg + 1]])
    np.savetxt(BASE_DIRECTORY + '/' + foldnameInversePower + 'xrange.txt',
               x_range,
               delimiter=',')
    np.savetxt(BASE_DIRECTORY + '/' + foldnameInversePower + 'yrange.txt',
               x_range,
               delimiter=',')

    print(x_range, 'x_range')
    print(y_range, 'y_range')

    pointset = []

    for x in x_range:
        for y in y_range:
            pointset.append((x % box_length, y % box_length))
    print(np.max(x_range % box_length))
    print(np.max(y_range % box_length))
    print(box_length)
    return pointset


def map_pointset_loop(foldname,
                      pointset,
                      quench,
                      coordarg,
                      use_minima_database=True,
                      minima_database_path=None):
    """ Checks a bunch of points if they match to a minimum by using a for loop
    """
    is_same_minimum_list = []
    resultlist = []
    ctol = 1e-1
    ndim = 2
    foldpath = BASE_DIRECTORY + '/' + foldname

    sysparams = load_params(foldpath)
    (hs_radii, initial_coords, box_length) = load_secondary_params(foldpath)

    minimum_coords = np.loadtxt(foldpath + '/coords_of_minimum.txt',
                                delimiter=',')
    minima_container = CheckSameMinimum(
        ctol,
        ndim,
        boxl=box_length,
        minimalist_max_len=2000,
        minima_database_location=minima_database_path,
        update_database=True)

    if use_minima_database == True:
        try:
            minima_container.minimalist = [
                minima_container.box_reshape_coords(x)
                for x in np.load(minima_database_path)
            ]
        except:
            print("warning no minima data found. generating")
            minima_container.minimalist = [
                minima_container.box_reshape_coords(minimum_coords)
            ]
    nfevlist = []
    nstepslist = []
    nhevlist = []
    for point in pointset:
        res = map_binary_inversepower(quench_steepest, foldname, point,
                                      coordarg)
        minima_container.add_minimum(res[0], point, res[2])
        nfevlist.append(res[3])
        nstepslist.append(res[4])
        nhevlist.append(res[5])

    print(np.average(nfevlist), 'number of function evaluations')
    print(np.average(nstepslist), 'number of steps')
    print(np.average(nstepslist), 'number of steps')
    print(np.average(nhevlist), "number of hessian evaluations")
    print(minima_container.orderparamlist)

    foldpathdata = foldpath + '/' + QUENCH_FOLDER_NAME
    os.makedirs(foldpathdata, exist_ok=True)

    minima_container.dump_map(foldpathdata)
    # print(minima_container.initial_coords_list)
    # print(minima_container.orderparamlist)
    # print(minima_container.orderparamlist)
    return is_same_minimum_list, resultlist


if __name__ == "__main__":
    foldnameInversePower = "ndim=2phi=0.9seed=0n_part=8r1=1.0r2=1.4rstd1=0.05rstd2=0.06999999999999999use_cell_lists=0power=2.5eps=1.0"
    coordarg = 0
    nmesh = 5
    pointset = construct_point_set_2d(foldnameInversePower, nmesh, 0.5,
                                      coordarg)
    # th = np.array(list(map(list, pointset))).T
    minima_database_path = BASE_DIRECTORY + '/' + foldnameInversePower + '/' + MINIMA_DATABASE_NAME
    res = map_pointset_loop(foldnameInversePower,
                            pointset,
                            quench_mixed_optimizer,
                            coordarg,
                            use_minima_database=True,
                            minima_database_path=minima_database_path)
    # print(res)
    # # boollist = np.loadtxt(BASE_DIRECTORY + '/' + foldnameInversePower + '/' + 'quench_results_fire.txt')
    # np.savetxt(BASE_DIRECTORY + '/' + foldnameInversePower + '/' + 'quench_results_mxopt.txt', boollist)
    # boollistreshaped = np.reshape(boollist, (nmesh, nmesh))
    # print(boollistreshaped)
    # print(boollist)
    # plt.imshow(boollistreshaped)
    # plt.show()
    # # print(reslist)
