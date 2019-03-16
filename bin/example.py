from actsims import simgen

"""
ACTSims example script
"""


# Use this version for s13-s15
#version = 'v4.0_mask_version_mr3c_20190215_pickupsub_190301'
version = 'v4.0_mask_version_mr3c_20190215_pickupsub_190303' # for s16

# our test data set
season, array, patch, freq = ('s13', 'pa1', 'deep1')
#season, array, patch, freq = ('s15', 'pa3', 'boss', 'f150')

# We initialize the sim generator with the mask version
simgen = simgen.SimGen(version=version)

# Signal + noise. Same note as above for stack of dichroic
# Only this function has been tested.
imap = simgen.get_sim(season, patch, array,sim_num=0)

# The output shape for get_sim is (nfreqs,nsplits,npol,Ny,Nx)
# The outputs will always have WCS that is "compatible" with the map-maker output, i.e. you can pad or slice as needed
# to what you are used to working with.
