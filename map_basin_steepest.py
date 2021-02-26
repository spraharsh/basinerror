"""
Maps the basin using the steepest descent algorithm. The accuracy of the mapping depends on the step size.
Note that our goal is to return the minimum given by the trajectory of the solution to d x/ d t = - \grad{E} (check the truth of this statement)
"""

from matplotlib.pyplot import sca
import numpy as np
import os


from pele.potentials.potential import potential
from quenches import quench_cvode_opt, quench_pycg_descent, quench_steepest
from params import load_params, load_secondary_params
from pele.potentials import InversePower
from params import BASE_DIRECTORY
import matplotlib.pyplot as plt
from pele.optimize._quench import modifiedfire_cpp, lbfgs_cpp
from quenches import quench_mixed_optimizer
from checksameminimum import CheckSameMinimum

# run data these parameter are stored separately so that they can be compared easily
# we don't want to delete any of these.
import optimizer_parameters as op
import optimizer_parameters_16 as op16
import optimizer_parameters_32 as op32
from pele.utils.cell_scale import get_ncellsx_scale



import write_run_data_to_file
from write_run_data_to_file import write_run_data_to_file


QUENCH_FOLDER_NAME = 'scratch'

# QUENCH_FOLDER_NAME = "cvode_high_tol_final"
# QUENCH_FOLDER_NAME = "fire_final"

# QUENCH_FOLDER_NAME = "lbfgs_m4_final"

# QUENCH_FOLDER_NAME = "lbfgs_m1_final"
# QUENCH_FOLDER_NAME = "CG_descent_final"
QUENCH_FOLDER_NAME = 'cvode_exact_initial'
# QUENCH_FOLDER_NAME = 'cvode_exact_lower'
MINIMA_DATABASE_NAME = 'minima_database_test_5.npy'

# things we need to define a run: foldername
# parameter dictionary
# name of the run

def map_binary_inversepower(foldername,
                            particle_coords,
                            optimizer, parameter_dict, random_coord_0=0, random_coord_1=-1):
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
    quench_coords = initial_coords.copy()
    print(initial_coords, 'initial')

    
    quench_coords[random_coord_0] = particle_coords[0]
    quench_coords[random_coord_1] = particle_coords[1]
    print(quench_coords, 'quench')
    # print(quench_coords, 'quench coords')
    # box length
    box_length = float(box_length)
    boxv = [box_length] * sysparams.ndim.value
    ncellx_scale = get_ncellsx_scale(hs_radii, boxv)
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
    try:
        ret = optimizer(quench_coords, potential, **parameter_dict)
    except:
        print(quench_coords, 'failed here')
        print(initial_coords, 'coords')
        print(len(quench_coords))
        raise Exception('failure')
    
    
    
    
    # ret = lbfgs_cpp(quench_coords, potential, tol=1e-8, M=1)

    # This exists because some runs don't have hessian evaluations
    try:
        ret['nhev']
    except:
        ret['nhev']=0
    coordarg = 0
    results = (ret.coords, ret.success, coordarg, ret.nfev, ret.nsteps, ret.nhev)
    # the reason the potential is being passed is because quench coords needs the potential to figure out what to do
    return results


def construct_point_set_2d(foldername, nmesh, boxlscale, random_coord_0=0, random_coord_1=1):
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
    # minimum_coords = np.loadtxt(foldpath + '/coords_of_minimum.txt',
    #                             delimiter=',')

    center = initial_coords

    # sets the length of the x and y range over which
    # the mesh is calculated
    xylength = boxlscale * box_length

    #  This division in cases is because linspace doesn't
    #  give the right answer for n = 1
    if (nmesh != 1):
        # x_range = np.linspace(minimum_coords[coordarg] - xylength/2,
        #   minimum_coords[coordarg] + xylength/2, nmesh)
        # y_range = np.linspace(minimum_coords[coordarg+1] - xylength/2,
        #   minimum_coords[coordarg+1] + xylength/2, nmesh)

        # position minimum in the center
        centered_x_of_min = center[random_coord_0] + box_length / 2
        centered_y_of_min = center[random_coord_1] + box_length / 2

        x_range = centered_x_of_min + np.linspace(0, box_length, nmesh)
        y_range = centered_y_of_min + np.linspace(0, box_length, nmesh)
    else:
        x_range = np.array([center[random_coord_0]])
        y_range = np.array([center[random_coord_1]])
    np.savetxt(BASE_DIRECTORY + '/' + foldnameInversePower + 'xrange.txt',
               x_range,
               delimiter=',')
    np.savetxt(BASE_DIRECTORY + '/' + foldnameInversePower + 'yrange.txt',
               x_range,
               delimiter=',')


    pointset = []

    for x in x_range:
        for y in y_range:
            pointset.append((x % box_length, y % box_length))
    return pointset



def map_pointset_loop(foldname,
                      pointset,
                      optimizer,
                      parameter_dict,
                      ctol=1e-2,
                      ndim=2,
                      use_minima_database=True,
                      minima_database_path=None, coord_arg_0=0, coord_arg_1=1):
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
    for index, point in enumerate(pointset):
        res = map_binary_inversepower(foldname, point,
                                      optimizer, parameter_dict, random_coord_0=coord_arg_0, random_coord_1=coord_arg_1)
        minima_container.add_minimum(res[0], point, res[2])
        # print(index)
        # print(minima_container.nrattlermin, 'nrattlermin')
        # print(minima_container.nfluidstates, 'nfluidstates')
        nfevlist.append(res[3])
        nstepslist.append(res[4])
        nhevlist.append(res[5])

    # print(np.average(nfevlist), 'number of function evaluations')
    # print(np.average(nstepslist), 'number of steps')
    # print(np.average(nstepslist), 'number of steps')
    # print(np.average(nhevlist), "number of hessian evaluations")
    
    # print(minima_container.orderparamlist)

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
    coord_arg_0 = -1
    coord_arg_1 = 1
    # set nmesh =1
    nmesh = 100
    pointset = construct_point_set_2d(foldnameInversePower, nmesh, 0.5, coord_arg_0, coord_arg_1)
    print(pointset)
    # th = np.array(list(map(list, pointset))).T

    # folders
    data_location = BASE_DIRECTORY + '/' + foldnameInversePower
    minima_database_path = data_location + '/' + MINIMA_DATABASE_NAME
    
    
    # defining parameter for run
    optimizer = quench_cvode_opt
    identification_tolerance = 1e-2
    parameter_dict = op.RUN_PARAMETERS_CVODE_EXACT_8
    optimizer_parameters = parameter_dict.copy()
    # important to remove name
    optimizer_parameters.pop('name', None)
    # run over the map
    res = map_pointset_loop(foldnameInversePower,
                            pointset,
                            optimizer,
                            optimizer_parameters,
                            ctol=identification_tolerance,
                            ndim=2,
                            use_minima_database=True,
                            minima_database_path=minima_database_path,
                            coord_arg_0=coord_arg_0, coord_arg_1=coord_arg_1)

    # print(res)
    # # boollist = np.loadtxt(BASE_DIRECTORY + '/' + foldnameInversePower + '/' + 'quench_results_fire.txt')
    # np.savetxt(BASE_DIRECTORY + '/' + foldnameInversePower + '/' + 'quench_results_mxopt.txt', boollist)
    # boollistreshaped = np.reshape(boollist, (nmesh, nmesh))
    # print(boollistreshaped)
    # print(boollist)
    # plt.imshow(boollistreshaped)
    # plt.show()
    # # print(reslist)


    sysparams = load_params(data_location)

    # save all parameters from run
    run_diagnostics = res[0]
    run_diagnostics['identification tolerance'] = identification_tolerance
    run_diagnostics['nmesh'] = nmesh
    run_diagnostics['ndim'] = sysparams.ndim.value
    run_diagnostics['nparticles'] = sysparams.n_part.value
    run_diagnostics['run data location'] = data_location + '/' + QUENCH_FOLDER_NAME

    

    opt_name = parameter_dict['name'].replace(" ", "_")
    # write run data to two different locations
    write_run_data_to_file(parameter_dict, run_diagnostics, folder_location=data_location, name=opt_name + '.yaml')
    write_run_data_to_file(parameter_dict, run_diagnostics,
                           folder_location='/home/praharsh/Dropbox/research/bv-libraries/basinerror/run_diagnostic_data',
                           name = opt_name + '.yaml')

