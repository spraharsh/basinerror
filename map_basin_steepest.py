"""
Maps the basin using the steepest descent algorithm. The accuracy of the mapping depends on the step size.
Note that our goal is to return the minimum given by the trajectory of the solution to d x/ d t = - \grad{E} (check the truth of this statement)
"""

from matplotlib.pyplot import sca
import numpy as np
import os

from pele.potentials.potential import potential
from quenches import quench_cvode_opt, quench_steepest
from params import load_params, load_secondary_params
from pele.potentials import InversePower
from params import BASE_DIRECTORY
import matplotlib.pyplot as plt
from pele.optimize._quench import modifiedfire_cpp, lbfgs_cpp
from quenches import quench_mixed_optimizer
from checksameminimum import CheckSameMinimum


# run data
import optimizer_parameters as op
import write_run_data_to_file
from write_run_data_to_file import write_run_data_to_file


QUENCH_FOLDER_NAME = 'fire2'
MINIMA_DATABASE_NAME = 'minima_database.npy'

# things we need to define a run: foldername
# parameter dictionary
# name of the run


def map_binary_inversepower(quench,
                            foldername,
                            particle_coords,
                            coordarg, optimizer, parameter_dict):
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
    #                              T=10,
    #                              step=1,
    #                              nsteps=100000,
    #                              conv_tol=1e-8,
    #                              tol=1e-6, rtol=1e-4, atol=1e-4)
    # ret = quench_steepest(
    #     potential,
    #     quench_coords,  # make sure right coords are being passed
    #     nsteps=2000000,
    #     stepsize=5e-3,  # for steepest descent step size should be small
    #     tol=1e-4)
    # ret = quench_cvode_opt(potential, quench_coords, tol=1e-6, rtol=1e-4, atol=1e-4)
    ret = optimizer(quench_coords, potential, **parameter_dict)
    
    
    # ret = lbfgs_cpp(quench_coords, potential, tol=1e-8, M=1)

    # This exists because some runs don't have hessian evaluations
    try:
        ret['nhev']
    except:
        ret['nhev']=0

    print(ret.nsteps, 'nsteps')
    results = (ret.coords, ret.success, coordarg, ret.nfev, ret.nsteps, ret.nhev)
    # the reason the potential is being passed is because quench coords needs the potential to figure out what to do
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
                      optimizer,
                      parameter_dict,
                      ctol=1e-2,
                      ndim=2,
                      use_minima_database=True,
                      minima_database_path=None):
    """ Checks a bunch of points if they match to a minimum by using a for loop
    """
    is_same_minimum_list = []
    resultlist = []
    foldpath = BASE_DIRECTORY + '/' + foldname

    sysparams = load_params(foldpath)
    (hs_radii, initial_coords, box_length) = load_secondary_params(foldpath)
    minimum_coords = np.loadtxt(foldpath + '/coords_of_minimum.txt',
                                delimiter=',')

    # Initialize CheckSameMinimum
    potential = InversePower(sysparams.power.value,
                             sysparams.eps.value,
                             use_cell_lists=False,
                             ndim=sysparams.ndim.value,
                             radii=hs_radii * 1.0,
                             boxvec=[box_length, box_length])
    minima_container = CheckSameMinimum(
        ctol,
        ndim,
        boxl=box_length,
        minimalist_max_len=2000,
        minima_database_location=minima_database_path,
        update_database=True, rattler_check=True, potential=potential, hs_radii=hs_radii)

    if use_minima_database == True:
        try:
            minima_container.minimalist = [
                minima_container.box_reshape_coords(x)
                for x in np.load(minima_database_path)
            ]
        except:
            print("warning no minima data found. generating")
            minima_container.minimalist = [
                # minima_container.box_reshape_coords(minimum_coords)
            ]
    nfevlist = []
    nstepslist = []
    nhevlist = []
    for point in pointset:
        res = map_binary_inversepower(quench_steepest, foldname, point,
                                      coordarg, optimizer, parameter_dict)
        minima_container.add_minimum(res[0], point, res[2])
        print(minima_container.nrattlermin, 'nrattlermin')
        print(minima_container.nfluidstates, 'nfluidstates')
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

    run_diagnostics = {}
    run_diagnostics['nfev'] = float(np.average(nfevlist))
    run_diagnostics['nhev'] = float(np.average(nhevlist))
    run_diagnostics['nsteps'] = float(np.average(nstepslist))
    
    
    # print(minima_container.initial_coords_list)
    # print(minima_container.orderparamlist)
    # print(minima_container.orderparamlist)
    return run_diagnostics, is_same_minimum_list, resultlist





if __name__ == "__main__":
    foldnameInversePower = "ndim=2phi=0.9seed=0n_part=8r1=1.0r2=1.4rstd1=0.05rstd2=0.06999999999999999use_cell_lists=0power=2.5eps=1.0"
    ndim = int()
    coordarg = 0
    nmesh = 10
    pointset = construct_point_set_2d(foldnameInversePower, nmesh, 0.5,
                                      coordarg)
    # th = np.array(list(map(list, pointset))).T


    # folders
    data_location = BASE_DIRECTORY + '/' + foldnameInversePower
    minima_database_path = data_location + '/' + MINIMA_DATABASE_NAME

    
    # defining parameter for run
    optimizer = modifiedfire_cpp
    identification_tolerance = 1e-2
    parameter_dict = op.RUN_PARAMETERS_MODIFIED_FIRE
    res = map_pointset_loop(foldnameInversePower,
                            pointset,
                            quench_mixed_optimizer,
                            coordarg,
                            optimizer,
                            parameter_dict,
                            ctol=identification_tolerance,
                            ndim=2,
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

    # save all parameters from run
    run_diagnostics = res[0]
    run_diagnostics['identification tolerance'] = identification_tolerance
    run_diagnostics['mesh size'] = identification_tolerance


    # write run data to two different locations
    write_run_data_to_file(parameter_dict, run_diagnostics, folder_location=data_location)
    write_run_data_to_file(parameter_dict, run_diagnostics, folder_location='/home/praharsh/Dropbox/research/bv-libraries/basinerror/run_diagnostic_data')
    