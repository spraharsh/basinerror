"""
Run parameters stored as dictionaries. These are going to be the base
on which runs shall be done. only supposed to be imported into
map_basin_steepest.py

Note that every time these are declared we want to write a separate set
to make sure that we always store ideal run parameters here
"""

RUN_PARAMETERS_MODIFIED_FIRE_32 = {
    "name": "modified_fire",
    "tol": 1e-8,
    "dtstart": 0.1,
    "dtmax": 1,
    "maxstep": 0.5,
    "Nmin": 5,
    "finc": 1.1,
    "fdec": 0.5,
    "fa": 0.99,
    "astart": 0.1
}


# note that CVODE wouldn't converge at higher tolerances at 1e-4 (need to check whether it's the rattler situation) 
RUN_PARAMETERS_CVODE_32 = {
    "name": "cvode",
    "tol": 1e-8,                   # tolerance with which minimum is identified
    "rtol": 1e-5,                # relative local tolerance of path
    "atol": 1e-5                 # relative absolute tolerance of path
}


# note that CVODE wouldn't converge at higher tolerances at 1e-4 (need to check whether it's the rattler situation) 
RUN_PARAMETERS_CVODE_32_TEST = {
    "name": "cvode_higher_tol_check",
    "tol": 1e-6,                   # tolerance with which minimum is identified
    "rtol": 1e-4,                # relative local tolerance of path
    "atol": 1e-4                 # relative absolute tolerance of path
}

# These parameters are for lower tolerance runs
# for figuring out exact basins
RUN_PARAMETERS_CVODE_EXACT_32 = {
    "name": "cvode_exact",
    "tol": 1e-8,                   # tolerance with which minimum is identified
    "rtol": 1e-10,                # relative local tolerance of path
    "atol": 1e-10                 # relative absolute tolerance of path
}

# run parameter at one order lower tolerance to figure out
RUN_PARAMETERS_CVODE_EXACT_LOWER_32 = RUN_PARAMETERS_CVODE_EXACT_32.copy()
RUN_PARAMETERS_CVODE_EXACT_LOWER_32['name'] = 'cvode_exact_lower'
RUN_PARAMETERS_CVODE_EXACT_LOWER_32['rtol'] = RUN_PARAMETERS_CVODE_EXACT_LOWER_32['rtol'] * 1e-1
RUN_PARAMETERS_CVODE_EXACT_LOWER_32['atol'] = RUN_PARAMETERS_CVODE_EXACT_LOWER_32['atol'] * 1e-1


RUN_PARAMETERS_CGDESCENT_32 = {
    "name": "CG_DESCENT",
    "tol": 1e-8,                   # tolerance with which minimum is identified
    "M":0
}
# we can just go with conv tol at this level and get good answers
# we can still be faster increasing it of course
# I also had to increase conv factor a bit because it messed up a few directions
RUN_PARAMETERS_MIXED_OPTIMIZER_32 = {
    # parameters CVODE
    "name": "mixed_optimizer_new",
    "tol": 1e-8,   # tolerance with which minimum is identified
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
    "conv_factor": 2.0
}

RUN_PARAMETERS_MIXED_OPTIMIZER_32_LOWER_TOL = {
    # parameters CVODE
    "name": "mixed_optimizer_new_lower_tol",
    "tol": 1e-8,   # tolerance with which minimum is identified
    "rtol": 1e-5,   # relative local tolerance of path
    "atol": 1e-5,   # relative absolute tolerance of path
    # parameters Newtons
    # number of steps at which a convergence step is run
    "T": 30,
    # Newton step size only make sense for stepsize 1
    "step": 1,
    # convexity tolerance i.e |\lambda_min/\lambda_max| <conv_tol
    "conv_tol":  1e-8,
    # convexity factor  addition to make it more convex
    "conv_factor": 2.0
}

RUN_PARAMETERS_MIXED_OPTIMIZER_32_LOWER_TOL_2 = {
    # parameters CVODE
    "name": "mixed_optimizer_new_lower_tol_2",
    "tol": 1e-8,   # tolerance with which minimum is identified
    "rtol": 1e-5,   # relative local tolerance of path
    "atol": 1e-5,   # relative absolute tolerance of path
    # parameters Newtons
    # number of steps at which a convergence step is run
    "T": 60,
    # Newton step size only make sense for stepsize 1
    "step": 1,
    # convexity tolerance i.e |\lambda_min/\lambda_max| <conv_tol
    "conv_tol":  1e-8,
    # convexity factor  addition to make it more convex
    "conv_factor": 2.0
}

RUN_PARAMETERS_CVODE_RTOL_1e_m6 = {
    "name": "cvode_exact_1e_m6",
    "tol": 1e-8,                   # tolerance with which minimum is identified
    "rtol": 1e-6,                # relative local tolerance of path
    "atol": 1e-6                 # relative absolute tolerance of path
}

RUN_PARAMETERS_CVODE_RTOL_1e_m7 = {
    "name": "cvode_exact_1e_m7",
    "tol": 1e-8,                   # tolerance with which minimum is identified
    "rtol": 1e-7,                # relative local tolerance of path
    "atol": 1e-7                 # relative absolute tolerance of path

}
RUN_PARAMETERS_MXOPT_RTOL_1e_m6 = {
    "name": "mxopt_1e_m6",
    "tol": 1e-8,                   # tolerance with which minimum is identified
    "rtol": 1e-6,                # relative local tolerance of path
    "atol": 1e-6,                 # relative absolute tolerance of path
    # parameters Newtons
    # number of steps at which a convergence step is run
    "T": 30,
    # Newton step size only make sense for stepsize 1
    "step": 1,
    # convexity tolerance i.e |\lambda_min/\lambda_max| <conv_tol
    "conv_tol":  1e-8,
    # convexity factor  addition to make it more convex
    "conv_factor": 2.0
}

RUN_PARAMETERS_MXOPT_RTOL_1e_m6 = {
    "name": "mxopt_1e_m6",
    "tol": 1e-8,                   # tolerance with which minimum is identified
    "rtol": 1e-6,                # relative local tolerance of path
    "atol": 1e-6,                 # relative absolute tolerance of path
    # parameters Newtons
    # number of steps at which a convergence step is run
    "T": 30,
    # Newton step size only make sense for stepsize 1
    "step": 1,
    # convexity tolerance i.e |\lambda_min/\lambda_max| <conv_tol
    "conv_tol":  1e-8,
    # convexity factor  addition to make it more convex
    "conv_factor": 2.0
}
RUN_PARAMETERS_MXOPT_RTOL_1e_m6_T100 = {
    "name": "mxopt_1e_m6_T100",
    "tol": 1e-8,                   # tolerance with which minimum is identified
    "rtol": 1e-6,                # relative local tolerance of path
    "atol": 1e-6,                 # relative absolute tolerance of path
    # parameters Newtons
    # number of steps at which a convergence step is run
    "T": 100,
    # Newton step size only make sense for stepsize 1
    "step": 1,
    # convexity tolerance i.e |\lambda_min/\lambda_max| <conv_tol
    "conv_tol":  1e-8,
    # convexity factor  addition to make it more convex
    "conv_factor": 2.0
}

RUN_PARAMETERS_MXOPT_RTOL_1e_m6_T300 = {
    "name": "mxopt_1e_m6_T300",
    "tol": 1e-8,                   # tolerance with which minimum is identified
    "rtol": 1e-6,                # relative local tolerance of path
    "atol": 1e-6,                 # relative absolute tolerance of path
    # parameters Newtons
    # number of steps at which a convergence step is run
    "T": 300,
    # Newton step size only make sense for stepsize 1
    "step": 1,
    # convexity tolerance i.e |\lambda_min/\lambda_max| <conv_tol
    "conv_tol":  1e-8,
    # convexity factor  addition to make it more convex
    "conv_factor": 2.0
}

RUN_PARAMETERS_MXOPT_RTOL_1e_m6_T200 = {
    "name": "mxopt_1e_m6_T300",
    "tol": 1e-8,                   # tolerance with which minimum is identified
    "rtol": 1e-6,                # relative local tolerance of path
    "atol": 1e-6,                 # relative absolute tolerance of path
    # parameters Newtons
    # number of steps at which a convergence step is run
    "T": 200,
    # Newton step size only make sense for stepsize 1
    "step": 1,
    # convexity tolerance i.e |\lambda_min/\lambda_max| <conv_tol
    "conv_tol":  1e-8,
    # convexity factor  addition to make it more convex
    "conv_factor": 2.0
}


RUN_PARAMETERS_MXOPT_RTOL_1e_m7 = {
    "name": "mxopt_1e_m7",
    "tol": 1e-8,                   # tolerance with which minimum is identified
    "rtol": 1e-7,                # relative local tolerance of path
    "atol": 1e-7,                 # relative absolute tolerance of path
    # parameters Newtons
    # number of steps at which a convergence step is run
    "T": 30,
    # Newton step size only make sense for stepsize 1
    "step": 1,
    # convexity tolerance i.e |\lambda_min/\lambda_max| <conv_tol
    "conv_tol":  1e-8,
    # convexity factor  addition to make it more convex
    "conv_factor": 2.0
}

# this happens to be LBFGS as written in pele.
# These parameters are NON STANDARD
# and not a good choice for generic optimization
# These parameters are the way they are because when they were
# chosen for generating more accurate basins
# of attraction fast in https://pubs.acs.org/doi/10.1021/jp312457a
RUN_PARAMETERS_LBFGS_M_4_32 = {
    "name": "LBFGS_M4",
    "tol": 1e-8,            # tolerance with which minimum is identified
    "M": 4,                    # size of stored memory
    # maximum size of step note that the step is made smaller
    "maxstep": 0.1,
    # uses sufficient increase condition (Check code for details)
    "maxErise": 1e-4,
    "H0": 0.1               # initial estimate of gradient
}

RUN_PARAMETERS_LBFGS_M_1_32 = {
    "name": "LBFGS_M1",
    "tol": 1e-8,            # tolerance with which minimum is identified
    "M": 1,                  # size of stored memory
    # maximum size of step note that the step is made smaller
    "maxstep": 0.1,
    # uses sufficient increase condition (Check code for details)
    "maxErise": 1e-4,
    "H0": 0.1               # initial estimate of gradient
}
