"""
Plots to do scaling analysis for error and function evaluations
This needs polishing
"""

import numpy as np
import matplotlib.pyplot as plt




# Old data correct basins

# n_vals = [8, 16, 32]
# error_cv_ht = [0.8, 1.6, 0.3]
# error_our_method = [2.7, 2.6, 3]
# error_fire = [11, 24, 24]
# error_CG = [15, 32, 30]
# error_LBFGS_M1 = [33, 55, 46]
# error_LBFGS_M4 = [24, 42, 36]

# nfev_cv_ht = [205, 352, 652]
# nfev_our_method = [153, 208, 268]
# nfev_fire = [198, 277, 491]
# nfev_CG = [88, 158, 312]
# nfev_LBFGS_M1 = [60, 98, 206]
# nfev_LBFGS_M4 = [52, 87, 167]

# New data random configs
n_vals = [8, 16, 32, 64]

# error_cv_ht = [0.6, 5.8, 8.8, 7.3]
error_our_method = [2.4, 5.8, 9.2, 13.3]
error_fire = [35, 67, 91, 99.4]
error_CG = [46, 76, 95, 99.7]
error_LBFGS_M1 =[81, 98, 100, 100]
error_LBFGS_M4 = [61, 88, 99, 100]

# nfev_cv_ht = [271, 507, 987, 1382]
nfev_our_method = [182, 276, 402, 835]
nfev_fire = [205, 318, 491, 878]
nfev_CG = [92, 179, 312, 600]
nfev_LBFGS_M1 = [72, 134, 206, 361]
nfev_LBFGS_M4 = [61, 112, 167, 283]






lw = 4
# plt.plot(n_vals, nfev_cv_ht, color="tab:blue", marker = 'o', label="CVODE high tol")
# Thickness change
plt.plot(n_vals, nfev_our_method, color="tab:orange", linewidth=lw, marker = 'o', label="our method")
plt.plot(n_vals, nfev_fire, color="tab:green", marker = 'o', label="fire")
plt.plot(n_vals, nfev_CG, color="tab:red", marker = 'o', label="CG")
plt.plot(n_vals, nfev_LBFGS_M1, color="tab:purple", marker = 'o', label="LBFGS_M1")
plt.plot(n_vals, nfev_LBFGS_M4, color="tab:brown", marker = 'o', label="LBFGS_M4")
plt.legend(loc="upper left")
plt.xlabel(r'Number of particles ($N$)')
plt.ylabel(r'gradient evaluations')
plt.savefig('nfevvsN_new_thick_transp.png', dpi=400, transparent=True)
plt.show()


# plt.plot(n_vals, error_cv_ht, color="tab:blue", marker = 'o', label="CVODE high tol")
# Thickness change
plt.plot(n_vals, error_our_method, color="tab:orange",  linewidth=lw,  marker = 'o', label="our method")
plt.plot(n_vals, error_fire, color="tab:green", marker = 'o', label="fire")
plt.plot(n_vals, error_CG, color="tab:red", marker = 'o', label="CG")
plt.plot(n_vals, error_LBFGS_M1, color="tab:purple", marker = 'o', label="LBFGS_M1")
plt.plot(n_vals, error_LBFGS_M4, color="tab:brown", marker = 'o', label="LBFGS_M4")
plt.legend(loc="right")
plt.xlabel(r'Number of particles ($N$)')
plt.ylabel(r'Error (%)')
plt.savefig('errorvsN_new_thick_transp.png', dpi=400, transparent=True)
plt.show()