import julia

from julia import Main
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
    box_length
    
    boxv = np.array([box_length]*sysparams.ndim.value)
    
    ncellx_scale = get_ncellsx_scale(radii, boxv)

    # potential = InversePower(sysparams.power.value,
    #                          sysparams.eps.value,
    #                          use_cell_lists=False,
    #                          ndim=sysparams.ndim.value,
    #                          radii=radii * 1.0,
    #                          boxvec=boxv)
    mxd = Main.Mixed_Descent(pele_wrapped_python_pot, Main.solver, Main.nls, coords, opt_param_dict['T'], opt_param_dict['rtol'], opt_param_dict['conv_tol'], opt_param_dict['tol'])

    

    pot = Main.pot.InversePower(
        sysparams.power.value,
        sysparams.eps.value,
        radii,
        ndim = sysparams.ndim.value,
        boxvec = boxv,
        use_cell_lists = False,
        ncellx_scale = cell_scale,
    )
    ppot = Main.PythonPotential(pot)
    

    try:
        Main.run_b(mxd, 10000)
    except:
        print("exception occured")
        # if exception occurs, treat as failure. this is for rattlers
        # not enough failures occur that it would make a difference to not just assume this never happens
        # but we shoudl switch this out
        # but jic
        return (quench_coords, False, 0, 0, 0, 0, 0)

    results = (mxd.optimizer.x0, mxd.converged, 0, mxd.iter_number, 0, 0, 0)
    return results


Main.include("mxopt_int.jl")
print("this works")
natoms = 64
radii_arr = Main.generate_radii(0, natoms, 1.0, 1.4, 0.05, 0.05 * 1.4)
dim = 2
phi = 0.9
power = 2.5
eps = 1

length_arr = Main.get_box_length(radii_arr, phi, dim)
coords = Main.generate_random_coordinates(length_arr, natoms, dim)


boxvec = np.array([length_arr, length_arr])

cell_scale = Main.utils.get_ncellsx_scale(radii_arr, boxvec)

th=Main.mxd.integrator.u

pele_wrapped_pot = Main.pot.InversePower(
    2.5,
    1.0,
    radii_arr,
    ndim = 2,
    boxvec = boxvec,
    use_cell_lists = False,
    ncellx_scale = cell_scale,
)



pele_wrapped_python_pot = Main.PythonPotential(pele_wrapped_pot)
mxd = Main.Mixed_Descent(pele_wrapped_python_pot, Main.solver, Main.nls, coords, 30, 10**(-5), 10^(-8)., 10**(-3))
Main.run_b(mxd, 2000)
print(mxd.integrator.u - coords)
print(pele_wrapped_python_pot)
print(th)
print("yay")





# nls = NewtonLinesearch(
#     lsolve_lsmr!,
#     p_energy,
#     p_gradient,
#     p_hessian,
#     backtracking_line_search!,
#     coords
# )







