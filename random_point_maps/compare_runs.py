"""
Script that compares dat between two runs
"""
from utils.cell_scale import get_box_length
from generate_points_random import get_hs_radii
from generate_points_random import SUB_FOLD_NAME
import yaml
import numpy as np
from params import BASE_DIRECTORY, load_params
# should rewrite this
from checksameminimum import CheckSameMinimum
from pele.potentials import InversePower
# faster loading and dumping with C backend
from yaml import CLoader, CDumper




def compare_runs_2d(fnames, foldpath, subfoldname, run_a, run_b, ctol):
    """ compares runs on points with names fnames between
        runs done with 2 minima find routines run_a and run_b

    Args:
       fnames filenames
       run_a run a 
       run_b run b

    Returns:
        percentage of minima that are the same

    """
    # TODO: replace this
    sysparams = load_params(foldpath)
    radii = get_hs_radii(foldpath, subfoldname)
    box_length = get_box_length(radii, sysparams.ndim.value, sysparams.phi.value)
    
    subfoldpath = foldpath + '/' + subfoldname
    data_path_a = subfoldpath + '/' + run_a
    data_path_b = subfoldpath + '/' + run_b
    potential = InversePower(sysparams.power.value,
                             sysparams.eps.value,
                             use_cell_lists=False,
                             ndim=sysparams.ndim.value,
                             radii=radii * 1.0,
                             boxvec=[box_length, box_length])
    minima_checker = CheckSameMinimum(ctol, dim=2, boxl=box_length,
                                      hs_radii=radii, potential=potential)

    
    same_minimum_check_l = []
    for fname in fnames:
        print(fname)
        minimum_a = np.loadtxt(data_path_a+ '/' + fname, delimiter=',')
        minimum_b = np.loadtxt(data_path_b+ '/' + fname, delimiter=',')
        boxed_minimum_a = minima_checker.box_reshape_coords(minimum_a)
        boxed_minimum_b = minima_checker.box_reshape_coords(minimum_b)
        # get first index that is not a rattler
        rattlers_exist, rattlers = minima_checker._find_rattlers(minimum_a)
        rattlers_dont_exist = np.count_nonzero(rattlers==0) == 0
        print(rattlers)
        # number of non rattlers to ensure we're not in a fluid state
        n_non_rattlers = np.count_nonzero(rattlers)
        # only do calculations if all particles aren't rattlers
        if n_non_rattlers != 0:
            first_non_rattler = (np.argwhere(rattlers!=0).T)[0,0]
            # TODO: rewrite the CheckSameMinimum function
            # load and make sure the particle being aligned is not a rattler later
            # we're choosing -1 because that's always going to be of radius 1.4
            #  particle
            aligned_minimum_b = minima_checker.align_structures(boxed_minimum_a, 
                                                            boxed_minimum_b, part_ind=first_non_rattler)
            same_minimum_check = minima_checker.check_same_structure(aligned_minimum_b,
                                                                 boxed_minimum_a, rattlers)
            same_minimum_check_l.append(same_minimum_check)

    fraction_same_minimum = np.mean(same_minimum_check_l)
    print(fraction_same_minimum)
    return fraction_same_minimum



def average_heuristics(foldpath, subfoldname, run_folder_name, fnames):
    subfoldpath = foldpath + '/' + subfoldname
    nfev_list = []
    nhev_list = []
    nsteps_list= []
    n_phase_1_list = []
    n_phase_2_list = []
    phase_heuristics_exist = True
    for fname in fnames:
        output_name = subfoldpath + '/' + run_folder_name + '/' + fname + 'heuristics.yaml'
        with open(output_name, 'r') as heur_file:
            heuristics = yaml.load(heur_file, Loader=CLoader)
            nfev_list.append(heuristics["nfev"])
            nhev_list.append(heuristics["nhev"])
            nsteps_list.append(heuristics["nsteps"])
            try:
                heuristics['n_phase_1']
                heuristics['n_phase_2']
            except:
                phase_heuristics_exist = False
            if phase_heuristics_exist:
                n_phase_1_list.append(heuristics['n_phase_1'])
                n_phase_2_list.append(heuristics['n_phase_2'])
            
    nfev_av = np.mean(nfev_list)
    nhev_av = np.mean(nhev_list)
    nsteps_av = np.mean(nsteps_list)
    if phase_heuristics_exist:
        n_phase_1_av = np.mean(n_phase_1_list)
        n_phase_2_av = np.mean(n_phase_2_list)
        av_dict = {'nfev_av':nfev_av, 'nhev_av':nhev_av, 'nsteps_av':nsteps_av, 'n_phase_1_av':n_phase_1_av,
                   'n_phase_2_av': n_phase_2_av}
    else:
        av_dict = {'nfev_av':nfev_av, 'nhev_av':nhev_av, 'nsteps_av':nsteps_av}
    # with open(subfoldpath + '/' + run_folder_name +  'average_heuristics.yaml', 'w') as heur_av_file:
    #     yaml.dump(av_dict, heur_av_file)
    return av_dict



def get_all_heuristics(fnames, foldpath, subfoldname, run_fold_name, correct_fold_name, ctol):
    """ Ouputs the percentage error and function evaluations for the optimimizer
    
    Args:
       foldpath ${1:arg1}
       subfoldpath ${2:arg2}
       run_fold_name ${3:arg3}
       fnames ${4:arg4}

    """
    fraction_wrong = 1 - (compare_runs_2d(fnames, foldpath, subfoldname, run_fold_name, correct_fold_name, ctol))
    percentage_wrong = fraction_wrong*100
    av_dict = average_heuristics(foldpath, subfoldname, run_fold_name, fnames)
    av_dict['percentage_wrong'] = float(percentage_wrong)
    # av_dict = {'percentage_wrong': float(percentage_wrong), 'nfev_av':float(nfev_av), 'nhev_av':float(nhev_av), 'nsteps_av':float(nsteps_av)}
    nfev_av = av_dict['nfev_av']
    nhev_av = av_dict['nhev_av']
    nsteps_av = av_dict['nsteps_av']
    try:
        n_phase_1_av = av_dict['n_phase_1_av']
        n_phase_2_av = av_dict['n_phase_2_av']
        print('n_phase_1' , n_phase_1_av)
        print('n_phase_2' , n_phase_2_av)
    except:
        print('no phase info')
    print("percentage wrong :" , fraction_wrong*100)
    print("nfev_av :", nfev_av)
    print("nhev_av :", nhev_av)
    print("nsteps_av :", nsteps_av)

    subfoldpath = foldpath + '/' + subfoldname
    with open(subfoldpath + '/' + run_fold_name +  '_heuristics.yaml', 'w') as heur_av_file:
        yaml.dump(av_dict, heur_av_file)
    
    



if __name__== "__main__":
    foldname = "ndim=2phi=0.9seed=0n_part=256r1=1.0r2=1.4rstd1=0.05rstd2=0.06999999999999999use_cell_lists=0power=2.5eps=1.0"
    sub_fold_name = SUB_FOLD_NAME
    fold_path = str(BASE_DIRECTORY+'/' + foldname)
    ensemble_size = int(500)
    fnames = list(map(str, range(ensemble_size)))
    # print(average_heuristics(sub_fold_path, opt_name, fnames))
    identification_tolerance = 1e-2    # comparison checks
    # opt_a = 'mixed_optimizer_new_lower_tol_2'
    opt_a = 'CVODE_high_tol_julia'
    # opt_b = 'cvode_exact_julia'
    #opt_a = 'cvode_julia'
    opt_a = 'mxopt_1e_m6_T100_julia_new'
    opt_b = 'cvode_exact_lower_julia'
    
    # print(compare_runs_2d(fnames, fold_path, sub_fold_name, opt_a, opt_b, identification_tolerance))
    get_all_heuristics(fnames, fold_path, sub_fold_name, opt_a, opt_b, identification_tolerance)