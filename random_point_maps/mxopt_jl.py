import julia

from julia import Main
Main.include("../quench_ode_julia/mxopt_int.jl")
import numpy as np

from utils.cell_scale import get_box_length, get_ncellsx_scale
from generate_points_random import SUB_FOLD_NAME, get_hs_radii
from params import BASE_DIRECTORY, load_params

def quench_single_mxopt_inverse_power_julia(coord_file_name, foldpath, sub_fold_name,
                                optimizer, opt_param_dict):
    """
    quenches a single system through mxopt
    Parameters
    ----------
    coord_file_name: string
        name of the path to the coordinates
    foldername: str
        folder definining the run
    sub_fold_name:
        name of subfolder where the run data is stored
    optimizer: optimizer
        quench
    opt_param_dict: dict
        dictionary of parameters for the optimizer
    """
    
    sysparams = load_params(foldpath)

    # path to quench coords
    quench_coords_path = foldpath + '/' + sub_fold_name + '/' + 'ensemble/' + coord_file_name
    quench_coords = np.loadtxt(quench_coords_path)
    radii = get_hs_radii(foldpath, sub_fold_name)

    box_length = get_box_length(radii, sysparams.ndim.value, sysparams.phi.value)
    
    boxv = np.array([box_length]*sysparams.ndim.value)
    
    ncellx_scale = get_ncellsx_scale(radii, boxv)

    # potential = InversePower(sysparams.power.value,
    #                          sysparams.eps.value,
    #                          use_cell_lists=False,
    #                          ndim=sysparams.ndim.value,
    #                          radii=radii * 1.0,
    #                          boxvec=boxv)
    pot = Main.pot.InversePower(
        sysparams.power.value,
        sysparams.eps.value,
        radii,
        ndim = sysparams.ndim.value,
        boxvec = boxv,
        use_cell_lists = False,
        ncellx_scale = ncellx_scale,
    )
    
    ppot = Main.PythonPotential(pot)
    nls = Main.NewtonLinesearch(ppot, quench_coords, opt_param_dict['tol'])
    mxd = Main.Mixed_Descent(ppot, Main.CVODE_BDF(), nls, quench_coords, opt_param_dict['T'], opt_param_dict['rtol'], opt_param_dict['conv_tol'], opt_param_dict['tol'])
    # Main.run_b(mxd, 10000)
    try:
        Main.run_b(mxd, 10000)
    except:
        print("exception occured")
        # if exception occurs, treat as failure. this is for rattlers
        # not enough failures occur that it would make a difference to not just assume this never happens
        # but we shoudl switch this out
        # but jic
        return (quench_coords, False, 0, 0, 0, 0, 0)
    results = (mxd.optimizer.x0, mxd.converged, mxd.n_g_evals, mxd.iter_number, mxd.n_h_evals, 0, 0)
    print(quench_coords_path)
    return results





# # Main.include("mxopt_int.jl")
# print("this works")
# natoms = 64
# radii_arr = Main.generate_radii(0, natoms, 1.0, 1.4, 0.05, 0.05 * 1.4)
# dim = 2
# phi = 0.9
# power = 2.5
# eps = 1

# length_arr = Main.get_box_length(radii_arr, phi, dim)
# coords = Main.generate_random_coordinates(length_arr, natoms, dim)


# boxvec = np.array([length_arr, length_arr])

# cell_scale = Main.utils.get_ncellsx_scale(radii_arr, boxvec)

# # th=Main.mxd.integrator.u

# pele_wrapped_pot = Main.pot.InversePower(
#     2.5,
#     1.0,
#     radii_arr,
#     ndim = 2,
#     boxvec = boxvec,
#     use_cell_lists = False,
#     ncellx_scale = cell_scale,
# )





# if __name__=="__main__":
#     solver = Main.CVODE_BDF()
#     pele_wrapped_python_pot = Main.PythonPotential(pele_wrapped_pot)
#     mxd = Main.Mixed_Descent(pele_wrapped_python_pot, solver, Main.nls, coords, 30, 10**(-5), 0., 10**(-3))