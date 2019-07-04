from actsims import simgen
from enlib import bench # for benchmarking

"""
ACTSims example script
"""


# You need to specify the mask version (for noise sims) first
# These are different for s16 and non-s16, sorry about that
version = 'v4.0_mask_version_mr3c_20190215_pickupsub_190301'
#version = 'v4.0_mask_version_mr3c_20190215_pickupsub_190303' # for s16

# our test data set
season, array, patch, freq = ('s13', 'pa1', 'deep1', 'f150')
#season, array, patch, freq = ('s15', 'pa3', 'boss', 'f150')

# We initialize the sim generator with the mask version

for theta in [0,15]:
    s09gen = simgen.Sehgal09Gen(version=version, eulers=(0,theta,0))

    # We can then get just signal = cmb + fg (and loop over season,patch,array,sim_num after the above initialization)
    with bench.show("signal"):
        s09gen.get_signal(season, patch, array, freq, sim_num=0)
    # Or get just cmb
    #s09gen.get_cmb(season, patch, array, freq, sim_num= 0)
    # Or get just foregrounds
    #s09gen.get_fg(season, patch, array, freq, sim_num=0)
    # Or get just noise. NOTE: You can't get this per frequency, but you get a stack of maps for the entire dichroic array
    with bench.show("noise"):
        s09gen.get_noise(season, patch, array,sim_num=0) #, freq, 0)
