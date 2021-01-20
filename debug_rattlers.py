"""
scratch space for debugging rattler configurations

TODO: write into tests later
"""


# coords correpsonding to issue
from pele.potentials import InversePower
from checksameminimum import CheckSameMinimum
from params import load_params, load_secondary_params
import numpy as np



debug_coords = [ 3.9324891118253564, 3.3781705856463096, 6.1520247054937478, 9.4293648218581847, 2.2382859620461653, 6.8602223824897655, 0.22373572115735679, 7.3563557948774783, 6.1043123397509182, 2.195355384152522, 2.0144892972106634, 8.8478063776409712, 3.7282664069461879, 5.4202113317897291, -0.67196297335207744, 5.4822416091893675, 6.0260496493211821, 4.6658280148182429, 7.1427659787530304, 7.2148532873658136, 1.4796105023022146, 4.594852954569471, 8.5626469148601583, 0.31125019065988851, 4.4776398646673909, 7.6827351326283706, 1.4531033229914481, 1.9008953836870657, 3.9645483931986578, 0.9724951066931532, 8.3663305823749905, 3.16224709493039 ]




ctol = 1e-2
ndim = 2
base_dir = '/home/praharsh/Dropbox/research/bv-libraries/basinerror/datainv'
foldnameInversePower = "ndim=2phi=0.9seed=0n_part=16r1=1.0r2=1.4rstd1=0.05rstd2=0.06999999999999999use_cell_lists=0power=2.5eps=1.0"
data_loc = base_dir + '/' + foldnameInversePower
minima_database_path = data_loc + 'minima_database.npy'
(hs_radii, initial_coords, box_length) = load_secondary_params(data_loc)
sysparams = load_params(data_loc)
potential = InversePower(sysparams.power.value,
             sysparams.eps.value,
             use_cell_lists=False,
             ndim=sysparams.ndim.value,
             radii=hs_radii * 1.0,
             boxvec=[box_length, box_length])
minima_container_debug = CheckSameMinimum(
        ctol,
        ndim,
        boxl=box_length,
        minimalist_max_len=2000,
        minima_database_location=minima_database_path,
        update_database=True, rattler_check=True, potential=potential, hs_radii=hs_radii)


# plot the rattlers





minima_container_debug._find_rattlers(np.array(debug_coords))

print('works')