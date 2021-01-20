"""
Code for plotting our lambdamin/lambdamax
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# this is obtained by writing the output of find_basin_minimum from CVODE
# 27f81661 to file and removing cruft from other output
data = np.loadtxt('output_n32.txt', delimiter=',')

print(data)
energies, conv_factors, lambda_min, lambda_max = data.T
print(data)
print(energies)
print(conv_factors)


fig, ax = plt.subplots()
ax.plot(-conv_factors, energies)
ax.set_xlabel(r'$\lambda_{min}/\lambda_{max}$')
ax.set_ylabel(r'Energy ($E$)')
plt.savefig('conv_vs_E.pdf')
plt.show()

fig, ax = plt.subplots()
ax.plot(lambda_min, energies, label=r'$\lambda_{min}$')
ax.plot(-lambda_max, energies, label=r'-$\lambda_{max}$')
left, bottom, width, height = [0.2, 0.2, 0.4, 0.4]
# ax2 = fig.add_axes([left, bottom, width, height])
# ax2.plot(-lambda_max, energies, color='green')
# ax2.set_xlabel(r'-$\lambda_{max}$')
# ax2.set_ylabel(r'Energy ($E$)')
# ax.plot(lambda_max, energies)
print(lambda_min)
ax.set_xlim(np.min(lambda_min), 0.1)  # decreasing time
# ax.set_xlabel(r'$\lambda_{min}$')
ax.set_xlabel(r'$\lambda_{min}$, $\lambda_{max}$')
ax.set_ylabel(r'Energy ($E$)')
plt.legend(loc=3)
plt.savefig('lmin_vs_E.pdf')
plt.show()


fig, ax = plt.subplots()
ax.plot(-lambda_min, label= r'$|\lambda_{min}|$')
ax.plot(conv_factors, label=r'$|\lambda_{min}/\lambda_{max}$|')
ax.plot(lambda_max, label=r'$|\lambda_{max}$|')
left, bottom, width, height = [0.4, 0.4, 0.4, 0.4]
# ax2 = fig.add_axes([left, bottom, width, height])
# ax2.plot(-lambda_min, label= r'$\lambda_{min}$')
# ax2.plot(conv_factors, label=r'$\lambda_{min}/\lambda_{max}$')
# ax2.plot(lambda_max, label=r'$\lambda_{max}$')
# ax2.set_xlabel(r'-$\lambda_{max}$')
# ax2.set_ylabel(r'Energy ($E$)')
# # ax2.set_yscale('log')
# ax2.set_xscale('symlog')
ax.set_xscale('symlog')
ax.set_xlabel(r'number of steps')
ax.set_ylabel(r'magnitude')
plt.legend()
plt.savefig('lamvsnsteps.pdf')
plt.show()