from actsims import simgen

version = 'v2.0_mask_version_mr3c_20190215_pickupsub_190301'
season, array, patch, freq = ('s15', 'pa3', 'deep56', 'f090')

simgen = simgen.SimGen(version=version)
simgen.get_signal(season, array, patch, freq, 0, 0)
simgen.get_cmb(season, array, patch, freq, 0, 0)
simgen.get_fg(season, array, patch, freq, 0, 0)
#simgen.get_noise(season, array, patch, freq, 0)

