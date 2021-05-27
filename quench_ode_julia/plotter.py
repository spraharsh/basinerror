import matplotlib.pyplot as plt
import numpy as np






tol = np.array([1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7])

nfev = np.array([261, 275, 358, 551, 916, 1501])
nhev = np.array([63, 74, 106, 175, 304, 508])
niter = np.array([206, 203, 252, 382, 629, 1029])
error = np.array([20, 19, 24, 27, 29, 31])




plt.xscale('log')
plt.yscale('log')

plt.xlabel('path tolerance (Diffeq.jl)')
plt.ylabel('evaluations/iterations')
plt.plot(tol, nfev, label= 'gradient evaluations')
plt.plot(tol, nhev, label= 'hessian evaluations')
plt.plot(tol, niter, label = 'iterations')
plt.plot(tol, error, label = 'percentage of points different from CVODE low tol (%)')
plt.title('gradient/hessian evaluations, iterations and error for QNDF')
plt.legend()
plt.savefig('accuracy_check_evaluations_QNBDF.pdf')
plt.show()
