from optimizer_parameters_32 import RUN_PARAMETERS_CGDESCENT_32, RUN_PARAMETERS_CVODE_32, RUN_PARAMETERS_CVODE_32_TEST, RUN_PARAMETERS_CVODE_EXACT_32, RUN_PARAMETERS_CVODE_EXACT_LOWER_32, RUN_PARAMETERS_CVODE_RTOL_1e_m6, RUN_PARAMETERS_CVODE_RTOL_1e_m7, RUN_PARAMETERS_LBFGS_M_1_32, RUN_PARAMETERS_LBFGS_M_4_32, RUN_PARAMETERS_MIXED_OPTIMIZER_32, RUN_PARAMETERS_MIXED_OPTIMIZER_32_LOWER_TOL, RUN_PARAMETERS_MIXED_OPTIMIZER_32_LOWER_TOL_2, RUN_PARAMETERS_MODIFIED_FIRE_32, RUN_PARAMETERS_MXOPT_RTOL_1e_m6, RUN_PARAMETERS_MXOPT_RTOL_1e_m6_T100, RUN_PARAMETERS_MXOPT_RTOL_1e_m7
from optimizer_parameters_16 import RUN_PARAMETERS_CGDESCENT_16, RUN_PARAMETERS_CVODE_16, RUN_PARAMETERS_CVODE_EXACT_16, RUN_PARAMETERS_CVODE_EXACT_LOWER_16, RUN_PARAMETERS_LBFGS_M_1_16, RUN_PARAMETERS_LBFGS_M_4_16, RUN_PARAMETERS_MIXED_OPTIMIZER_T_30_16, RUN_PARAMETERS_MODIFIED_FIRE_16
from optimizer_parameters import *
from utils.cell_scale import get_box_length, get_ncellsx_scale
from generate_points_random import SUB_FOLD_NAME, get_hs_radii
from pele.optimize._quench import modifiedfire_cpp, lbfgs_cpp
import numpy as np
import os
import yaml
from pele.potentials import InversePower
from params import BASE_DIRECTORY, load_params
from quenches import quench_LBFGS, quench_cvode_opt, quench_mixed_optimizer, quench_pycg_descent, quench_steepest
from pele.potentials import Harmonic

from timeit import default_timer as timer
def quench_single_inverse_power(coord_file_name, foldpath, sub_fold_name,
                                optimizer, opt_param_dict):
    """
    figures out the minimum correspoding to a set of particle coords
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
    
    boxv = [box_length]*sysparams.ndim.value
    ncellx_scale = get_ncellsx_scale(radii, boxv)

    potential = InversePower(sysparams.power.value,
                             sysparams.eps.value,
                             use_cell_lists=False,
                             ndim=sysparams.ndim.value,
                             radii=radii * 1.0,
                             boxvec=boxv)

    # quench_coords = [1., 1., 1.]
    # potential = Harmonic(np.array([0.,0., 0.]), 1)
    

    try:
        ret = optimizer(quench_coords, potential, **opt_param_dict)
    except:
        print("exception occured")
        # if exception occurs, treat as failure. this is for rattlers
        # not enough failures occur that it would make a difference to not just assume this never happens
        # but we shoudl switch this out
        # but jic
        return (quench_coords, False, 0, 0, 0, 0, 0)

    # This exists because some runs don't have hessian evaluations
    try:
        ret['nhev']
    except:
        ret['nhev']=0

    # mixed optimizer statistics
    try:
        ret['n_phase_1']
    except:
        ret['n_phase_1'] = 0
    # mixed optimizer statistics
    try:
        ret['n_phase_2']
    except:
        ret['n_phase_2'] = 0

    results = (ret.coords, ret.success, ret.nfev, ret.nsteps, ret.nhev, ret.n_phase_1, ret.n_phase_2)
    print(quench_coords_path)
    return results



def quench_multiple(foldpath, sub_fold_name, fnames, output_dir,
                    optimizer, opt_param_dict):
    """
    does multiple quenches and saves results for each quench.
    result coordinates are saved as output_dir/fname for each fname
    and run heuristics as output_dir/fname_heuristics.yml
    Parameters
    ----------
    foldpath: name of configuration folder
        folder as defined for parameters
    sub_fold_name: folder name
        name of folder containing data
    output_dir: folder name
       folder containing output data
    fnames: List(str)
        list of paths in ensemble that need to be changed
    optimizer: optimizer
        quench function
    opt_param_dict: dict
        dictionary of optimizer parameters
    """
    total_quenches = len(fnames)
    os.makedirs(foldpath + '/' + sub_fold_name + '/' + output_dir, exist_ok=True)
    for index, fname in enumerate(fnames):
        print('quenching ' + str(index) + ' out of ' + str(total_quenches))
        
        results = quench_single_inverse_power(fname, foldpath, sub_fold_name, optimizer, opt_param_dict)
        final_coords, success, nfev, nsteps, nhev, n_phase_1, n_phase_2 = results
        heuristics_dict = {'success': success, 'nfev': nfev, 'nsteps':nsteps, 'nhev':nhev, 'n_phase_1':n_phase_1, 'n_phase_2':n_phase_2}
        print(heuristics_dict)
        heuristics_dict.update(opt_param_dict)
        # save accordingly
        output_name = foldpath + '/' + sub_fold_name + '/' + output_dir + '/' + fname
        
        print(output_name)
        print("is output ")
        np.savetxt(output_name, final_coords)
        with open(output_name + 'heuristics.yaml', 'w') as heur_file:
            yaml.dump(heuristics_dict, heur_file)





if __name__== "__main__":
    foldname = "ndim=2phi=0.9seed=0n_part=8r1=1.0r2=1.4rstd1=0.05rstd2=0.06999999999999999use_cell_lists=0power=2.5eps=1.0"
    foldpath = str(BASE_DIRECTORY+'/' + foldname)
    ensemble_size = int(1e4)
    # quench = quench_cvode_opt
    # opt_params = RUN_PARAMETERS_CVODE_EXACT_LOWER_32
    # opt_name= opt_params['name']
    # opt_params.pop('name', None)
    # fnames = list(map(str, range(ensemble_size)))
    # quench_multiple(foldpath, SUB_FOLD_NAME, fnames,opt_name, quench, opt_params)
    # print(fnames)
    # quench = quench_cvode_opt
    opt_params = RUN_PARAMETERS_CVODE_32
    # opt_name= opt_params['name']
    # opt_params.pop('name', None)
    # fnames = list(map(str, range(ensemble_size)))
    # quench_multiple(foldpath, SUB_FOLD_NAME, fnames,opt_name, quench, opt_params)
    # print(fnames)
    # quench=quench_mixed_optimizer
    # opt_params = RUN_PARAMETERS_MXOPT_RTOL_1e_m6
    # opt_name= opt_params['name']
    # opt_params.pop('name', None)
    # fnames = list(map(str, range(ensemble_size)))
    # quench_multiple(foldpath, SUB_FOLD_NAME, fnames,opt_name, quench, opt_params)
    quench=quench_cvode_opt
    # opt_params = RUN_PARAMETERS_CVODE_32
    opt_name= opt_params['name']
    opt_params.pop('name', None)
    fnames = list(map(str, range(ensemble_size)))
    fnames = fnames[:1]
    start = timer()
    quench_multiple(foldpath, SUB_FOLD_NAME, fnames,opt_name, quench, opt_params)
    end = timer()
    print(end-start)
