"""
Run parameters stored as dictionaries. These are going to be the base
on which runs shall be done. only supposed to be imported into
map_basin_steepest.py

Note that every time these are declared we want to write a separate set
to make sure that we always store ideal run parameters here
"""

RUN_PARAMETERS_MODIFIED_FIRE_8 = {
    "name": "modified fire",
    "tol": 1e-6,
    "dtstart": 0.1,
    "dtmax": 1,
    "maxstep": 0.5,
    "Nmin": 5,
    "finc": 1.1,
    "fdec": 0.5,
    "fa": 0.99,
    "astart": 0.1
}

RUN_PARAMETERS_CVODE_8 = {
    "name": "CVODE_high_tol",
    "tol": 1e-6,                   # tolerance with which minimum is identified
    "rtol": 1e-4,                # relative local tolerance of path
    "atol": 1e-4                 # relative absolute tolerance of path
}

RUN_PARAMETERS_JULIA_8 = {
    "name": "julia_high_tol",
    "tol": 1e-6,                   # tolerance with which minimum is identified
    "rtol": 1e-1,                # relative local tolerance of path
    "atol": 1e-1                 # relative absolute tolerance of path
}


# These parameters are for lower tolerance runs
# for figuring out exact basins
RUN_PARAMETERS_CVODE_EXACT_8 = {
    "name": "CVODE_exact",
    "tol": 1e-6,                   # tolerance with which minimum is identified
    "rtol": 1e-8,                # relative local tolerance of path
    "atol": 1e-8                 # relative absolute tolerance of path
}

# run parameter at one order lower tolerance to figure out
RUN_PARAMETERS_CVODE_EXACT_LOWER_8 = RUN_PARAMETERS_CVODE_EXACT_8.copy()
RUN_PARAMETERS_CVODE_EXACT_LOWER_8['name'] = 'CVODE_exact_lower_tol'
RUN_PARAMETERS_CVODE_EXACT_LOWER_8['rtol'] = RUN_PARAMETERS_CVODE_EXACT_LOWER_8['rtol'] * 1e-1
RUN_PARAMETERS_CVODE_EXACT_LOWER_8['atol'] = RUN_PARAMETERS_CVODE_EXACT_LOWER_8['atol'] * 1e-1


RUN_PARAMETERS_CGDESCENT_8 = {
    "name": "CG_DESCENT",
    "tol": 1e-6,                   # tolerance with which minimum is identified
    "M":0
}



RUN_PARAMETERS_MIXED_OPTIMIZER_8 = {
    # parameters CVODE
    "name": "mixed optimizer",
    "tol": 1e-6,   # tolerance with which minimum is identified
    "rtol": 1e-4,   # relative local tolerance of path
    "atol": 1e-4,   # relative absolute tolerance of path
    # parameters Newtons
    # number of steps at which a convergence step is run
    "T": 10,
    # Newton step size only make sense for stepsize 1
    "step": 1,
    # convexity tolerance i.e |\lambda_min/\lambda_max| <conv_tol
    "conv_tol":  1e-2,
    # convexity factor  addition to make it more convex
    "conv_factor": 2,
}


# this happens to be LBFGS as written in pele.
# These parameters are NON STANDARD
# and not a good choice for generic optimization
# These parameters are the way they are because when they were
# chosen for generating more accurate basins
# of attraction fast in https://pubs.acs.org/doi/10.1021/jp312457a
RUN_PARAMETERS_LBFGS_M_4_8 = {
    "name": "LBFGS_M4",
    "tol": 1e-6,            # tolerance with which minimum is identified
    "M": 4,                    # size of stored memory
    # maximum size of step note that the step is made smaller
    "maxstep": 0.1,
    # uses sufficient increase condition (Check code for details)
    "maxErise": 1e-4,
    "H0": 0.1               # initial estimate of gradient
}

RUN_PARAMETERS_LBFGS_M_1_8 = {
    "name": "LBFGS_M1",
    "tol": 1e-6,            # tolerance with which minimum is identified
    "M": 1,                  # size of stored memory
    # maximum size of step note that the step is made smaller
    "maxstep": 0.1,
    # uses sufficient increase condition (Check code for details)
    "maxErise": 1e-4,
    "H0": 0.1               # initial estimate of gradient
}


RUN_PARAMETERS_MIXED_OPTIMIZER_T_30_8 = {
    # parameters CVODE
    "name": "mixed_optimizer_new",
    "tol": 1e-6,   # tolerance with which minimum is identified
    "rtol": 1e-4,   # relative local tolerance of path
    "atol": 1e-4,   # relative absolute tolerance of path
    # parameters Newtons
    # number of steps at which a convergence step is run
    "T": 30,
    # Newton step size only make sense for stepsize 1
    "step": 1,
    # convexity tolerance i.e |\lambda_min/\lambda_max| <conv_tol
    "conv_tol":  1e-8,
    # convexity factor  addition to make it more convex
    "conv_factor": 2,
}