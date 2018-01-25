# lensSims

This is a simple set of tools to return sims of noise and signal for CMB studies with ACTPol.

The sims can be accessed using simTools.getActpolSim(iterationNum, simType, patch)

where you should set:
* iterationNum is an integer from 0 to 127 
* simType = 'cmb' or 'noise' ('foregrounds' to come)
* patch is a string set to, for instance (see the inputParams/CoaddAll*dict file for the current full list),
    patch = 'deep1_s13_pa1_f150'
* Alternatively, you can pass the info separately, e.g:
    season = 's13', 
    pa = 'pa1', 
    freqGHz = 150,
    region = 'deep5'
 
See the bin/simpleTest.py script for a simple invocation.


To install: pip install -e . --user

Current issues / to-do:
* Noise sims and signal sims currently are getting slightly different sizes (originating because one is obtained using Flipper, the other with enlib).  There is a mechanism for dealing with this (the coordsEpsilonArcmin flag) but the values for each patch are not yet determined.
* Boss-N is to come.
* simType = 'foregrounds' to be implemented
* For lensing, we need a second set of cmb sims lensed by the same phi field; these have not been made yet, so cmbSim = 0 has been set as the default internally in the code for now.
* Currently only sims of noise for full coadds are generated (not splits)
* FFTs for noise sims currently take a long time, as they are using numpy.fft rather than fftw.  This can be replaced (or better, replace with enlib versions)
