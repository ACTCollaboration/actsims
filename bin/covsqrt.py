from orphics import maps,io,stats
from enlib import bench,array_ops
from pixell import enmap,utils
import numpy, numpy as np
import os
from tilec import covtools


# mtypes = [
#     # 's15_deep56_pa3_unsmoothed_pscov.fits',
#     's15_deep56_pa3_smoothed_pscov.fits']

mtypes = [
    's13_deep6_pa1_unsmoothed_pscov.fits',
    's13_deep6_pa1_smoothed_pscov.fits']

for mtype in mtypes:
    icov = enmap.read_map(mtype)
    cov = enmap.enmap(maps.symmat_from_data(icov).to_array(),icov.wcs)
    
    print(mtype, " multipow")
    covsqrt = enmap.enmap(enmap.multi_pow(cov.copy(),0.5),icov.wcs)
    enmap.write_map(mtype.replace('.fits','_multipow_covsqrt.fits'),covsqrt)


    print(mtype, " arrayops")
    covsqrt = array_ops.eigpow(cov.copy(), 0.5, axes = [0,1])
    enmap.write_map(mtype.replace('.fits','_arrayops_covsqrt.fits'),covsqrt)


    print(mtype, " multipow_noiqu")
    ncomps = cov.shape[0]
    for i in range(ncomps):
        for j in range(ncomps):
            if i==j: continue
            if (i==3 and j==0) or (i==0 and j==3): continue
            cov[i,j] = 0
    covsqrt = enmap.enmap(enmap.multi_pow(cov.copy(),0.5),icov.wcs)
    enmap.write_map(mtype.replace('.fits','_multipow_noiqu_covsqrt.fits'),covsqrt)
    
    print(mtype, " arrayops_noiqu")
    covsqrt = array_ops.eigpow(cov.copy(), 0.5, axes = [0,1])
    enmap.write_map(mtype.replace('.fits','_arrayops_noiqu_covsqrt.fits'),covsqrt)
