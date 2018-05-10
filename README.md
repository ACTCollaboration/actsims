# actsims

This is a simple set of tools to return sims of noise and signal for CMB studies with ACTPol.


## To run
The sims can be accessed using `simTools.getActpolSim(iterationNum, simType, psa)`

where you should set:
* `iterationNum` is an integer from 0 to 2047 
* `simType = 'cmb'`,  `'noise'`, or `'foregrounds'`
* `psa` stands for "patch, season, array", and is a string set to, for instance,  `patch = 'deep1_s13_pa1'` (see the `inputParams/templateInput.dict` file for the current full list)
* Alternatively, you can pass the info separately, e.g:
    `season = 's13', 
    pa = 'pa1', 
    region = 'deep5'`
    
    See the `bin/simpleTest.py` script for a simple invocation.

The code will return an `enmap` object of stacked maps; the shape is `[N_freq, 3, Ny, Nx]`, where the first index is 1 or 2 depending on whether the array is dichroic, and the second index is for `T, Q, U`


## Installation
To install: `pip install -e . --user`

Make sure that you have pulled a version of enlib after Jan 25 2018 (Sigurd merged in Mat's fix for an [issue](https://github.com/amaurea/enlib/issues/34) on that date)

### Feynman installation

In the `actsims` directory, do the following:

```
ln -s /mnt/act3/users/msyriac/sim_data/beams beams
ln -s /mnt/act3/users/msyriac/sim_data/signal/v0.2 cmb_sim_data
ln -s /mnt/act3/users/msyriac/sim_data/noise/v0.2/ templateData
```


## Changes in April 26 edition
* New power spectrum templates for all the maps, that give reasonable power spectra, see http://www.cita.utoronto.ca/~engelen/plots/sim_vs_nullmaps.pdf
* If using a dichroic array, sims at both frequencies are returned; the correlations between the frequencies are included in  the foregrounds, and the same is true in principle for the noise.  
* Using enmap everywhere (no flipper), all code is now self-contained (no code from cmblens)
* Pixel window functions now being applied
* Signal sims have the modulation & aberration applied, and can be labelled separately for CMB and phi maps 

## Current issues / to-do:
* Diagonlize the big (6, 6, Ny, Nx) matrix for the full covariance, and then test the sims coming out of this pipeline for that template.
* Add a test routine to check whether the code returns the same realization (across platfroms and python module versions).
* CMB sims are awaiting a new run, to fix an issue with the aberration effect and with the random seeds.
* Need to do lensing nulls for these.
