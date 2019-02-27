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
with bench.show("imports"):
    from tilec import covtools
import argparse

# Parse command line
parser = argparse.ArgumentParser(description='Make covsqrt, generate some test sims, make verification plots.')
parser.add_argument("version", type=str,help='A prefix for a unique version name')
parser.add_argument("model", type=str,help='Name of a datamodel specified in soapack.interfaces.')
parser.add_argument("--mask-version", type=str,  default="180323",help='Mask version')
parser.add_argument("--mask-kind", type=str,  default="binary_apod",help='Mask kind')
parser.add_argument("--mask-patch", type=str,  default=None,help='Mask patch')
parser.add_argument("--mask-pad", type=int,  default=None,
                    help='Mask additional padding. No padding is applied to the extracted mask if any specified.')
parser.add_argument("--extract-mask", type=str,  default=None,
                    help='Make sims on the big mask but do all the analysis on an extract of this version.')
parser.add_argument("--covsqrt-kind", type=str,default="arrayops",help='Method for covsqrt.')
parser.add_argument("--binary-percentile", type=float,  default=10.,help='Binary percentile for sim masking.')
parser.add_argument("--season", type=str,help='Season')
parser.add_argument("--array", type=str,help='Array')
parser.add_argument("--patch", type=str,help='Patch')
parser.add_argument("-n", "--nsims",     type=int,  default=10,help="Number of sims.")
parser.add_argument("-r", "--radial-fit-annulus",     type=int,  default=20,help="Bin width for azimuthal averaging.")
parser.add_argument("-d", "--dfact",     type=int,  default=8,help="Downsample factor.")
parser.add_argument("-a", "--aminusc", action='store_true',help='Whether to use the auto minus cross estimator.')
parser.add_argument("--no-write", action='store_true',help='Do not write any FITS to disk.')
parser.add_argument("--no-off", action='store_true',help='Null the off-diagonals.')
parser.add_argument("--no-prewhiten", action='store_true',help='Do not prewhiten spectra before smoothing. Use this flag for Planck.')
parser.add_argument("--overwrite", action='store_true',help='Overwrite an existing version.')
parser.add_argument("--debug", action='store_true',help='Debug plots.')
args = parser.parse_args()
coadd = not(args.aminusc)
nsims = args.nsims
if args.mask_patch is None: mask_patch = args.patch
else: mask_patch = args.mask_patch
if args.binary_percentile < 1e-3: bp = None
else: bp = args.binary_percentile
if args.dfact == 0: 
    smooth = False
else: 
    smooth = True
    dfact = (args.dfact,args.dfact)

# Make version tag
version = args.version
other_keys={'mask_version':args.mask_version}
for key in other_keys.keys():
    version += ("_"+key+"_"+str(other_keys[key]))
#####


# Get file name convention
pout,cout,sout = noise.get_save_paths(args.model,version,coadd,
                                      season=args.season,patch=args.patch,array=args.array,
                                      mkdir=True,overwrite=args.overwrite)
# Get data model
mask = sints.get_act_mr3_crosslinked_mask(mask_patch,
                                          version=args.mask_version,
                                          kind=args.mask_kind,
                                          season=args.season,array=args.array+"_f150",
                                          pad=args.mask_pad)
with bench.show("data model"):
    dm = sints.models[args.model](region=mask)

# Get a NoiseGen model
if args.extract_mask is not None:
    emask = sints.get_act_mr3_crosslinked_mask(mask_patch,version=args.extract_mask,kind=args.mask_kind,season=args.season,array=args.array+"_f150")
    eshape,ewcs = emask.shape,emask.wcs
else:
    emask = mask
ngen = noise.NoiseGen(version=version,model=args.model,extract_region=emask,ncache=0)

# Get arrays from array

with bench.show("load data"):
    splits = dm.get_splits(season=args.season,patch=args.patch,arrays=dm.array_freqs[args.array],srcfree=True)
    if args.debug: noise.plot(pout+"_splits",splits)
    ivars = dm.get_splits_ivar(season=args.season,patch=args.patch,arrays=dm.array_freqs[args.array])
modlmap = splits.modlmap()
with bench.show("n2d"):
    n2d_flat = noise.get_n2d_data(splits,ivars,mask,coadd_estimator=coadd,flattened=True,plot_fname=pout+"_n2d_flat" if args.debug else None)
del splits
radial_pairs = [(0,0),(1,1),(2,2),(3,3),(4,4),(5,5),(0,3),(3,0)] if not(args.no_prewhiten) else []
if smooth:
    n2d_flat_smoothed = noise.smooth_ps(n2d_flat.copy(),dfact=dfact,
                                        radial_pairs=radial_pairs,
                                        plot_fname=pout+"_n2d_flat_smoothed" if args.debug else None,
                                        radial_fit_annulus = args.radial_fit_annulus)
else:
    n2d_flat_smoothed = n2d_flat.copy()
del n2d_flat
if args.no_off: n2d_flat_smoothed = noise.null_off_diagonals(n2d_flat_smoothed)

with bench.show("covsqrt"):
    covsqrt = noise.get_covsqrt(n2d_flat_smoothed,args.covsqrt_kind)
del n2d_flat_smoothed
ngen.save_covsqrt(covsqrt,season=args.season,patch=args.patch,array=args.array,coadd=coadd)

if nsims>0:
    bin_edges = np.arange(40,8000,40)
    p1ds = []
    for i in range(nsims):
        with bench.show("print"):
            print("Sim %d of %d ..." % (i+1,nsims))
        with bench.show("simgen"):
            sims = ngen.generate_sim(season=args.season,patch=args.patch,array=args.array,seed=i,binary_percentile=bp)
            print(sims.nbytes/1024./1024./1024., " GB", sims.shape, sims.dtype)
        if args.extract_mask is not None: 
            ivars2 = enmap.extract(ivars,eshape,ewcs)
            modlmap = enmap.modlmap(eshape,ewcs)
        else:
            ivars2 = ivars

        if args.debug and i==0: noise.plot(pout+"_sims",sims)
        enmap.write_map("%s_trial_sim_seed_%d.fits" % (sout,i) ,sims)
        n2d_sim = noise.get_n2d_data(sims,ivars2,emask,coadd_estimator=coadd,flattened=False,plot_fname=pout+"_n2d_sim" if args.debug else None)
        del sims
        cents,op1ds_sim = noise.get_p1ds(n2d_sim,modlmap,bin_edges)
        p1ds.append(op1ds_sim.copy().reshape(-1))
    p1dstats = stats.get_stats(np.array(p1ds))

    del covsqrt

    # For verification
    splits = dm.get_splits(season=args.season,patch=args.patch,arrays=dm.array_freqs[args.array],srcfree=True)

    if args.extract_mask is not None: 
        splits = enmap.extract(splits,eshape,ewcs)

    n2d_data = noise.get_n2d_data(splits,ivars2,emask,coadd_estimator=coadd)
    cents,p1ds_data = noise.get_p1ds(n2d_data,modlmap,bin_edges)
    corr = noise.corrcoef(n2d_data)
    del n2d_data
    bin_edges = np.arange(40,10000,100)
    cents2,c1ds_data = noise.get_p1ds(corr,modlmap,bin_edges)
    noise.plot_corrcoeff(cents2,c1ds_data,plot_fname=pout)

    nfreqs = len(dm.array_freqs[args.array])
    noise.compare_ps(cents,p1dstats['mean'].reshape((nfreqs*3,nfreqs*3,cents.size)),p1ds_data,plot_fname="%s_compare" % (pout),err=p1dstats['errmean'].reshape((nfreqs*3,nfreqs*3,cents.size)))
