from actsims import simgen
from enlib import bench # for benchmarking
"""
ACTSims example script
"""
# You need to specify the mask version (for noise sims) first
# These are different for s16 and non-s16, sorry about that
#version = 'v4.0_mask_version_mr3c_20190215_pickupsub_190301'
version = 'v4.0_mask_version_mr3c_20190215_pickupsub_190303' # for s16
# our test data set
season, array, patch, freq = ('s16', 'pa2', 'cmb', 'f150')
#season, array, patch, freq = ('s13', 'pa1', 'deep1', 'f150')
#season, array, patch, freq = ('s15', 'pa3', 'boss', 'f150')
# We initialize the sim generator with the mask version
simgen = simgen.SimGen(version=version)
# We can then get just signal = cmb + fg (and loop over season,patch,array,sim_num after the above initialization)
# with bench.show("signal"):
#     simgen.get_signal(season, patch, array, freq, sim_num=0,mask_patch="patch006")
# Or get just cmb
#imap = simgen.get_cmb(season, patch, array, freq, sim_num= 0,mask_patch="patch006")
# Or get just foregrounds
simgen.get_fg(season, patch, array, freq, sim_num=0,mask_patch="patch006")
# Or get just noise. NOTE: You can't get this per frequency, but you get a stack of maps for the entire dichroic array
#with bench.show("noise"):
#    simgen.get_noise(season, patch, array,sim_num=0) #, freq, 0)
# Or signal + noise. Same note as above for stack of dichroic
# @Steve: hopefully this is the only function you will need
#imap = simgen.get_sim(season, patch, array,sim_num=0,mask_patch="patch006")
