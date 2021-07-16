from numpy.core.records import array
from optimizer_parameters_32 import RUN_PARAMETERS_CGDESCENT_32, RUN_PARAMETERS_CVODE_32, RUN_PARAMETERS_CVODE_32_TEST, RUN_PARAMETERS_CVODE_EXACT_32, RUN_PARAMETERS_CVODE_EXACT_LOWER_32, RUN_PARAMETERS_CVODE_RTOL_1e_m6, RUN_PARAMETERS_CVODE_RTOL_1e_m7, RUN_PARAMETERS_LBFGS_M_1_32, RUN_PARAMETERS_LBFGS_M_4_32, RUN_PARAMETERS_MIXED_OPTIMIZER_32, RUN_PARAMETERS_MIXED_OPTIMIZER_32_LOWER_TOL, RUN_PARAMETERS_MIXED_OPTIMIZER_32_LOWER_TOL_2, RUN_PARAMETERS_MODIFIED_FIRE_32, RUN_PARAMETERS_MXOPT_RTOL_1e_m6, RUN_PARAMETERS_MXOPT_RTOL_1e_m6_T100, RUN_PARAMETERS_MXOPT_RTOL_1e_m7
from optimizer_parameters_16 import RUN_PARAMETERS_CGDESCENT_16, RUN_PARAMETERS_CVODE_16, RUN_PARAMETERS_CVODE_EXACT_16, RUN_PARAMETERS_CVODE_EXACT_LOWER_16, RUN_PARAMETERS_LBFGS_M_1_16, RUN_PARAMETERS_LBFGS_M_4_16, RUN_PARAMETERS_MIXED_OPTIMIZER_16, RUN_PARAMETERS_MIXED_OPTIMIZER_T_30_16, RUN_PARAMETERS_MODIFIED_FIRE_16
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
# from quench_ode_julia.interface import ode_julia_naive
# from mxopt_jl import quench_single_mxopt_inverse_power_julia
# from assigner_jl import quench_single_CVODE_inverse_power_julia
from timeit import default_timer as timer
from multiprocessing import Pool

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

    quench_coords = np.array([2.352688223090913, 6.3280683094189225, 3.8090762744749256, 3.5065782011179385, 1.6908106110661583, 5.917323888564995, 1.0861062799713277, 4.261870546129833, 3.7445472114103735, 5.645740209262758, 6.292846747674822, 4.991692333803862, 3.0515165149950967, 0.4055276571722094, 2.2961213822804694, 4.991117299492139])
    radii = np.array([1.0339553713017888, 1.0414206741450018, 0.9823496299849702, 0.9932573064034739, 1.4410631952243176, 1.420813509559459, 1.404546328398384, 1.3923687830442797])

    quench_coords = [9.337647221542596, 1.5456510115859738, 9.90887183280318, 10.359109125537435, 0.4565747264079197, 6.260760163366896, 11.615697890514323, 12.293774616628575, 10.348709383183378, 1.6059793439488075, 1.4761039517598815, 1.0276707949222579, 10.032942247470187, 1.354080713122294, 10.826106545008862, 2.3783636187898707, 4.032236132596851, 2.5371569145787833, 11.284765973194554, 8.460162809464745, 7.579049049185771, 8.173050435342091, 9.494651056572234, 7.083387877400749, 6.160286900406435, 7.781370134554123, 10.222472937606318, 0.6345486477165924, 6.236228510824097, 8.593603020694431, 10.880113121189947, 2.1844985603481306, 12.127990212130813, 10.124036287647558, 6.533896105572216, 8.03690012368885, 3.9632426322036567, 10.287283606649998, 6.534793951918366, 5.064563240852603, 6.192243612316144, 0.8812982127394248, 6.199886221361277, 6.050671029744465, 5.844756686997425, 8.321805765435386, 8.64558240534, 1.224504799376647, 12.257515762225177, 10.013743200726655, 2.2917380712380284, 1.0337186826593456, 8.81245791227988, 12.005307391106015, 12.5705787449805, 8.419549901249782, 7.70589767982475, 6.191008990916649, 0.8168852579970988, 6.07961010280596, 4.940458887243859, 12.614904306548436, 2.791610189390051, 11.150533589443071]
    radii = [1.0339553713017888, 1.0414206741450018, 0.9823496299849702, 0.9932573064034739, 1.0293308537316554, 1.0148667925424708, 1.003247377427417, 0.9945491307459141, 0.9742894804583339, 1.0787165101068494, 0.9655546436087151, 0.9618598091794771, 1.019874120460908, 1.0405814811253438, 0.9826822697860601, 0.990621365902742, 1.2874920573110558, 1.226344508853804, 1.5593362982949166, 1.4153785427279777, 1.3918003728486945, 1.357912248407023, 1.4799593507258924, 1.3937968577259667, 1.4195626376462929, 1.4077995174664062, 1.3749481539730772, 1.4331600025504532, 1.4210163750316778, 1.3466126058232804, 1.499613400646502, 1.4285870771873195]

    radii = np.array(radii)
    quench_coords = np.array(quench_coords)
    
    boxv = [box_length]*sysparams.ndim.value
    boxv = np.array([6.502221171875915, 6.502221171875915])
    boxv = [12.917825730504669, 12.917825730504669]
    boxv = np.array(boxv)
    ncellx_scale = get_ncellsx_scale(radii, boxv)

    potential = InversePower(sysparams.power.value,
                             sysparams.eps.value,
                             use_cell_lists=False,
                             ndim=sysparams.ndim.value,
                             radii=radii * 1.0,
                             boxvec=boxv)
    
    boxv = [box_length]*sysparams.ndim.value
    # ncellx_scale = get_ncellsx_scale(radii, boxv)

    print(potential.getEnergy(quench_coords))
    ret = optimizer(quench_coords, potential, **opt_param_dict)
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

    print(ret.coords-quench_coords)
    print(opt_param_dict)

    results = (ret.coords, ret.success, ret.nfev, ret.nsteps, ret.nhev, ret.n_phase_1, ret.n_phase_2)
    print(quench_coords_path)
    print(results)
    return results



def quench_multiple(foldpath, sub_fold_name, fnames, output_dir,
                    optimizer, opt_param_dict, quench = quench_single_inverse_power):
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
        
        results = quench(fname, foldpath, sub_fold_name, optimizer, opt_param_dict)
        final_coords, success, nfev, nsteps, nhev, n_phase_1, n_phase_2 = results
        heuristics_dict = {'success': success, 'nfev': nfev, 'nsteps':nsteps, 'nhev':nhev, 'n_phase_1':n_phase_1, 'n_phase_2':n_phase_2}
        heuristics_dict.update(opt_param_dict)
        # save accordingly
        output_name = foldpath + '/' + sub_fold_name + '/' + output_dir + '/' + fname
        np.savetxt(output_name, final_coords)
        with open(output_name + 'heuristics.yaml', 'w') as heur_file:
            yaml.dump(heuristics_dict, heur_file)

def quench_multiple_parallel(foldpath, sub_fold_name, fnames, output_dir,
                             optimizer, opt_param_dict, quench = quench_single_inverse_power, nprocs = 1):
    """
    does multiple quenches and saves results for each quench.
    result coordinates are saved as output_dir/fname for each fname
    and run heuristics as output_dir/fname_heuristics.yml. This function 
    runs the jobs in parallel. This  does not work with julia quenches
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
    def single_run(fname):
        results = quench(fname, foldpath, sub_fold_name, optimizer, opt_param_dict)
        final_coords, success, nfev, nsteps, nhev, n_phase_1, n_phase_2 = results
        heuristics_dict = {'success': success, 'nfev': nfev, 'nsteps':nsteps, 'nhev':nhev, 'n_phase_1':n_phase_1, 'n_phase_2':n_phase_2}
        heuristics_dict.update(opt_param_dict)
        # save accordingly
        output_name = foldpath + '/' + sub_fold_name + '/' + output_dir + '/' + fname
        np.savetxt(output_name, final_coords)
        with open(output_name + 'heuristics.yaml', 'w') as heur_file:
            yaml.dump(heuristics_dict, heur_file)
    os.makedirs(foldpath + '/' + sub_fold_name + '/' + output_dir, exist_ok=True)
    with Pool(nprocs) as pool:
        pool.map(single_run, fnames)







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
    opt_params = RUN_PARAMETERS_MIXED_OPTIMIZER_32
    opt_params['name'] += '_test'

    print(opt_params)
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
    # quench= quench_mixed_optimizer
    print(opt_params)
    quench = quench_mixed_optimizer
    opt_params = RUN_PARAMETERS_MXOPT_RTOL_1e_m6_T100
    opt_name= opt_params['name']
    opt_params.pop('name', None)
    opt_params['nsteps'] = 3000
    # opt_params['sparse'] = True
    fnames = list(map(str, range(ensemble_size)))
    
    fnames = fnames[:1]
    
    start = timer()
    #quench_multiple(foldpath, SUB_FOLD_NAME, fnames,opt_name, quench, opt_params, quench=quench_single_mxopt_inverse_power_julia)
    quench_multiple(foldpath, SUB_FOLD_NAME, fnames,opt_name, quench, opt_params, quench=quench_single_inverse_power)
    end = timer()
    print("time: " + str(end - start))