"""
This contains additional tools for making sure two minima are actually different.

1. Check what the gradient and hessian are at the two minima
2. check Newtons from 1 point to the other 
"""


from map_basin_steepest import MINIMA_DATABASE_NAME
from params import BASE_DIRECTORY
import numpy as np
# bunch of useful functions here for periodic boundary conditions
from checksameminimum import CheckSameMinimum
from params import load_params, load_secondary_params
from pele.potentials import InversePower


# change according to which minima you want to compare
m1_arg = 148
m2_arg = 159
m3_arg = 194
m4_arg = 195



foldnameInversePower = "ndim=2phi=0.9seed=0n_part=16r1=1.0r2=1.4rstd1=0.05rstd2=0.06999999999999999use_cell_lists=0power=2.5eps=1.0"
minima_database_path = BASE_DIRECTORY + '/' + foldnameInversePower + '/' + MINIMA_DATABASE_NAME
th = np.load(minima_database_path)


foldpath = BASE_DIRECTORY + '/' + foldnameInversePower
sysparams = load_params(foldpath)
(hs_radii, initial_coords, box_length) = load_secondary_params(foldpath)
ctol = 1e-3
ndim = 2
potential = InversePower(sysparams.power.value,
                             sysparams.eps.value,
                             use_cell_lists=False,
                             ndim=sysparams.ndim.value,
                             radii=hs_radii * 1.0,
                             boxvec=[box_length, box_length])
minima_container = CheckSameMinimum(
    ctol,
    ndim,
    boxl=box_length,
    minimalist_max_len=2000,
    minima_database_location=minima_database_path,
    update_database=True, rattler_check=True, potential=potential, hs_radii=hs_radii)
print(box_length, 'box_length')



box_length = float(box_length)
minimum_coords = np.loadtxt(foldpath + '/coords_of_minimum.txt',
                            delimiter=',')
m1 = minima_container.box_reshape_coords(th[m1_arg])
m2 = minima_container.box_reshape_coords(th[m2_arg])
m3 = minima_container.box_reshape_coords(th[m3_arg])
m4 = minima_container.box_reshape_coords(th[m4_arg])
print(minima_container.boxl)


# print(th)
# print(m1)
# print(m2)
# print(m3)
# print(m4)
print(m3)
print(m4)
dv = m3 - m4
print(dv)
print(np.where(dv>box_length/2, box_length-dv, dv))

# dist
# print(m1-m2)
# print(m3-m4)
# a = m3-m4
# print(a[1,1]-box_length, 'weird point')
# print(m3-m4)
