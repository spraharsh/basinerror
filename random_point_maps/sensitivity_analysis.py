"""
Script that compares dat between two runs
"""
from generate_points_random import SUB_FOLD_NAME
import yaml
import numpy as np
from params import BASE_DIRECTORY




def compare_runs(fnames, run_a, run_b):
    """ compares runs on points with names fnames between
        runs done with 2 minima find routines run_a and run_b

    Args:
       fnames filenames
       run_a run a 
       run_b run b

    Returns:
        percentage of minima that are the same

    """
    return 0



def average_heuristics(run_folder, fnames):
    nfev_list = []
    nhev_list = []
    nsteps_list= []
    for fname in fnames:
        output_name = run_folder + '/' + fname + 'heuristics.yaml'
        with open(output_name, 'r') as heur_file:
            heuristics = yaml.load(heur_file, Loader=yaml.FullLoader)
            nfev_list.append(heuristics["nfev"])
            nhev_list.append(heuristics["nhev"])
            nsteps_list.append(heuristics["nsteps"])
            
    nfev_av = np.mean(nfev_list)
    nhev_av = np.mean(nhev_list)
    nsteps_av = np.mean(nsteps_list)
    av_dict = {'nfev_av':nfev_av, 'nhev_av':nhev_av, 'nsteps_av':nsteps_av}
    with open(run_folder + '/' + 'average_heuristics.yaml', 'w') as heur_av_file:
        yaml.dump(av_dict, heur_av_file)
    return nfev_av, nhev_av, nsteps_av



if __name__== "__main__":
    foldname = "ndim=2phi=0.9seed=0n_part=8r1=1.0r2=1.4rstd1=0.05rstd2=0.06999999999999999use_cell_lists=0power=2.5eps=1.0"
    opt_name = 'cvode_exact'
    sub_fold_name = SUB_FOLD_NAME
    foldpath = str(BASE_DIRECTORY+'/' + foldname + '/' + sub_fold_name + '/' + opt_name)
    ensemble_size = int(5e3)
    fnames = list(map(str, range(ensemble_size)))
    print(average_heuristics(foldpath, fnames))
    print(foldpath)
    




