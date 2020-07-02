from __future__ import print_function
from orphics import maps,io,cosmology,stats
from pixell import enmap,lensing,curvedsky as cs,utils
import numpy as np
import os,sys
from enlib import bench

shape,wcs = enmap.fullsky_geometry(res=1.0 * utils.arcmin)
lmax = 8000

theory = cosmology.default_theory()
ps = cosmology.enmap_power_from_orphics_theory(theory,lensed=False)
with bench.show('lens'):
    uTquMap, lTquMap, pMap = lensing.rand_map((3,)+shape, wcs, ps, lmax=lmax, output="ulp", verbose=True,
                                              phi_seed = 0,
                                              seed = 0, dtype=np.float32,delta_theta = 60 * utils.degree)


"""
dtheta = 10 deg : 273 s , 21 GB
dtheta = 20 deg : 271 s , 22 GB
dtheta = 40 deg : 268 s , 26 GB
dtheta = 60 deg : 267 s , 28 GB

OMP_NUM_THREADS=10
dtheta = 10 deg : 
dtheta = 20 deg : 364 s , 13 GB
dtheta = 40 deg : 
dtheta = 60 deg : 357 s , 20 GB

OMP_NUM_THREADS=20
dtheta = 10 deg : 
dtheta = 20 deg : 300 s , 15 GB
dtheta = 40 deg : 
dtheta = 60 deg : 

"""
