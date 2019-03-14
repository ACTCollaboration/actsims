from actsims import simgen

"""

cmbSeedInd = 0
fgSeedInd = 1
phiSeedInd = 2
noiseSeedInd = 3

phiSeed = (0, 0, 2, i)
cmbSeed = (set, 0, 0, i)
fgseed = (set, 0, 1, i)

"""

version = 'v1.0_mask_version_mr3c_20190215_pickupsub_190301'
#season, array, patch, freq = ('s15', 'pa2', 'deep56', 'f090')
season, array, patch, freq = ('s13', 'pa1', 'deep1', 'f150')

simgen = simgen.SimGen(version=version)
# simgen.get_signal(season, array, patch, freq, 0, 0)
# simgen.get_cmb(season, array, patch, freq, 0, 0)
# simgen.get_fg(season, array, patch, freq, 0, 0)
# simgen.get_noise(season, patch, array,sim_num=0) #, freq, 0)
imap = simgen.get_sim(season, patch, array,sim_num=0)
from orphics import io
io.hplot(imap,"simtest")
print(imap.shape)

