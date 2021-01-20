"""
Plots to do scaling analysis for error and function evaluations
This needs polishing
"""


n_vals = [8, 16, 32]

error_cv_ht = [0.8, 1.6, 0.3]
error_our_method = [2.7, 2.6, 3]
error_fire = [11, 24, 24]
error_CG = [15, 32, 30]
error_LBFGS_M1 = [33, 55, 46]
error_LBFGS_M4 = [24, 42, 36]


nfev_cv_ht = [205, 352, 652]
nfev_our_method = [153, 208, 268]
nfev_fire = [198, 277, 491]
nfev_CG = [88, 158, 312]
nfev_LBFGS_M1 = [60, 98, 206]
nfev_LBFGS_M4 = [52, 87, 167]

plt.plot(n_vals, nfev_cv_ht, color="tab:blue", marker = 'o', label="CVODE high tol")
plt.plot(n_vals, nfev_our_method, color="tab:orange", marker = 'o', label="our method")
plt.plot(n_vals, nfev_fire, color="tab:green", marker = 'o', label="fire")
plt.plot(n_vals, nfev_CG, color="tab:red", marker = 'o', label="CG")
plt.plot(n_vals, nfev_LBFGS_M1, color="tab:purple", marker = 'o', label="LBFGS_M1")
plt.plot(n_vals, nfev_LBFGS_M4, color="tab:brown", marker = 'o', label="LBFGS_M4")
plt.legend(loc="upper left")
plt.xlabel(r'Number of particles ($N$)')
plt.ylabel(r'gradient evaluations')
plt.savefig('nfevvsN.pdf')
plt.show()


plt.plot(n_vals, error_cv_ht, color="tab:blue", marker = 'o', label="CVODE high tol")
plt.plot(n_vals, error_our_method, color="tab:orange", marker = 'o', label="our method")
plt.plot(n_vals, error_fire, color="tab:green", marker = 'o', label="fire")
plt.plot(n_vals, error_CG, color="tab:red", marker = 'o', label="CG")
plt.plot(n_vals, error_LBFGS_M1, color="tab:purple", marker = 'o', label="LBFGS_M1")
plt.plot(n_vals, error_LBFGS_M4, color="tab:brown", marker = 'o', label="LBFGS_M4")
plt.legend(loc="upper left")
plt.xlabel(r'Number of particles ($N$)')
plt.ylabel(r'Error (%)')
plt.savefig('errorvsN.pdf')
plt.show()