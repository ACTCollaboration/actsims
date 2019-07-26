"""
This script can be used to make a covsqrt and a few trial sims.
"""
from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
from pixell import enmap,enplot
import numpy as np
import os,sys
from actsims import noise,utils
from soapack import interfaces as sints
from enlib import bench
from orphics import io,stats
import matplotlib.pyplot as plt
from tilec import covtools
import argparse

# Parse command line
parser = argparse.ArgumentParser(description='Make covsqrt, generate some test sims, make verification plots.')
parser.add_argument("--mask-pad", type=int,  default=None,
                    help='Mask additional padding. No padding is applied to the extracted mask if any specified.')
parser.add_argument("--covsqrt-kind", type=str,default="arrayops",help='Method for covsqrt.')
parser.add_argument("--array", type=str,help='Array')
parser.add_argument("--mask-patch", type=str,help='Patch')
parser.add_argument("-n", "--nsims",     type=int,  default=3,help="Number of sims.")
parser.add_argument("-r", "--radial-fit-annulus",     type=int,  default=20,help="Bin width for azimuthal averaging.")
parser.add_argument("--no-off", action='store_true',help='Null the off-diagonals.')
parser.add_argument("-d", "--dfact",     type=int,  default=8,help="Downsample factor.")
args = parser.parse_args()
coadd = True
nsims = args.nsims
if args.dfact == 0: 
    smooth = False
else: 
    smooth = True
    dfact = (args.dfact,args.dfact)

# Make version tag
version = "debug_s16"
other_keys={}
for key in other_keys.keys():
    version += ("_"+key+"_"+str(other_keys[key]))
#####
season = "s16"
mask_patch = args.mask_patch
patch = "cmb"
mask_version = "mr3c_20190215_pickupsub_190303"
model = "act_mr3"

# Get file name convention
pout,cout,sout = noise.get_save_paths(model,version,coadd,
                                      season=season,patch=patch,array=args.array,
                                      mkdir=True,overwrite=True,mask_patch=mask_patch)
# Get data model
mask = sints.get_act_mr3_crosslinked_mask(mask_patch,
                                          version=mask_version,
                                          kind="binary_apod",
                                          season=season,array=args.array+"_f150",
                                          pad=args.mask_pad)
# print(mask.shape)
mask = mask[400:-400,5000:-5000] 
mask *= enmap.apod(mask,400)
# print(mask.shape)
noise.plot(pout+"_mask",mask,grid=True)
dm = sints.models[model](region=mask)

# Get a NoiseGen model
emask = mask
ngen = noise.NoiseGen(version=version,model=model,extract_region=emask,ncache=1,verbose=True)

# Get arrays from array

splits = dm.get_splits(season=season,patch=patch,arrays=dm.array_freqs[args.array],srcfree=True)
ivars = dm.get_splits_ivar(season=season,patch=patch,arrays=dm.array_freqs[args.array])
noise.plot(pout+"_splits",splits)
noise.plot(pout+"_ivars",ivars)

modlmap = splits.modlmap()
n2d_flat = noise.get_n2d_data(splits,ivars,mask,coadd_estimator=coadd,flattened=True,plot_fname=pout+"_n2d_flat")
del splits

radial_pairs = [(0,0),(1,1),(2,2),(3,3),(4,4),(5,5),(0,3),(3,0)] 
if smooth:
    n2d_flat_smoothed = noise.smooth_ps(n2d_flat.copy(),dfact=dfact,
                                        radial_pairs=radial_pairs,
                                        plot_fname=pout+"_n2d_flat_smoothed",
                                        radial_fit_annulus = args.radial_fit_annulus)
else:
    n2d_flat_smoothed = n2d_flat.copy()
if args.no_off: n2d_flat_smoothed = noise.null_off_diagonals(n2d_flat_smoothed)
del n2d_flat

covsqrt = noise.get_covsqrt(n2d_flat_smoothed,args.covsqrt_kind)
del n2d_flat_smoothed
ngen.save_covsqrt(covsqrt,season=season,patch=patch,array=args.array,coadd=coadd,mask_patch=mask_patch)

if nsims>0:
    bin_edges = np.arange(40,8000,40)
    p1ds = []
    for i in range(nsims):
        print("Sim %d of %d ..." % (i+1,nsims))
        with bench.show("simgen"):
            sims = ngen.generate_sim(season=season,patch=patch,array=args.array,seed=i,binary_percentile=None,mask_patch=mask_patch)
            print(sims.nbytes/1024./1024./1024., " GB", sims.shape, sims.dtype)
        ivars2 = ivars

        noise.plot(pout+"_sims",sims)
        n2d_sim = noise.get_n2d_data(sims,ivars2,emask,coadd_estimator=coadd,flattened=False,plot_fname=pout+"_n2d_sim" if (i==0) else None)
        del sims
        cents,op1ds_sim = noise.get_p1ds(n2d_sim,modlmap,bin_edges)
        p1ds.append(op1ds_sim.copy().reshape(-1))
    p1dstats = stats.get_stats(np.array(p1ds))

    del covsqrt

    # For verification
    splits = dm.get_splits(season=season,patch=patch,arrays=dm.array_freqs[args.array],srcfree=True)

    n2d_data = noise.get_n2d_data(splits,ivars2,emask,coadd_estimator=coadd)
    cents,p1ds_data = noise.get_p1ds(n2d_data,modlmap,bin_edges)
    corr = noise.corrcoef(n2d_data)
    del n2d_data
    bin_edges = np.arange(40,10000,100)
    cents2,c1ds_data = noise.get_p1ds(corr,modlmap,bin_edges)
    noise.plot_corrcoeff(cents2,c1ds_data,plot_fname=pout)

    nfreqs = len(dm.array_freqs[args.array])
    noise.compare_ps(cents,p1dstats['mean'].reshape((nfreqs*3,nfreqs*3,cents.size)),p1ds_data,plot_fname="%s_compare" % (pout),err=p1dstats['errmean'].reshape((nfreqs*3,nfreqs*3,cents.size)))
