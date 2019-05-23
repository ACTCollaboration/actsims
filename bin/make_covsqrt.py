"""
This script can be used to make a covsqrt and a few trial sims.
"""
from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
from pixell import enmap,enplot,fft
import numpy as np
import os,sys
from actsims import noise
from soapack import interfaces as sints
from enlib import bench
from orphics import io,stats
import matplotlib.pyplot as plt
from tilec import covtools
import argparse

# Parse command line
parser = argparse.ArgumentParser(description='Make covsqrt, generate some test sims, make verification plots.')
parser.add_argument("version", type=str,help='A prefix for a unique version name')
parser.add_argument("model", type=str,help='Name of a datamodel specified in soapack.interfaces.')
parser.add_argument("--do-only-filter-noise", action='store_true',help='Do not do noise sim templates. Instead just do unflattened filter noise.')
parser.add_argument("--mask-version", type=str,  default="padded_v1",help='Mask version')
parser.add_argument("--mask-kind", type=str,  default="binary_apod",help='Mask kind')
parser.add_argument("--mask-patch", type=str,  default=None,help='Mask patch')
parser.add_argument("--mask-pad", type=int,  default=None,
                    help='Mask additional padding. No padding is applied to the extracted mask if any specified.')
parser.add_argument("--extract-mask", type=str,  default=None,
                    help='Make sims on the big mask but do all the analysis on an extract of this version.')
parser.add_argument("--covsqrt-kind", type=str,default="arrayops",help='Method for covsqrt.')
parser.add_argument("--season", type=str,help='Season')
parser.add_argument("--array", type=str,help='Array')
parser.add_argument("--patch", type=str,help='Patch')
parser.add_argument("--rlmin",     type=int,  default=300,help="Minimum ell.")
parser.add_argument("-n", "--nsims",     type=int,  default=10,help="Number of sims.")
parser.add_argument("-r", "--radial-fit-annulus",     type=int,  default=20,help="Bin width for azimuthal averaging.")
parser.add_argument("-d", "--dfact",     type=int,  default=8,help="Downsample factor.")
parser.add_argument("-a", "--aminusc", action='store_true',help='Whether to use the auto minus cross estimator.')
parser.add_argument("--no-write", action='store_true',help='Do not write any FITS to disk.')
parser.add_argument("--calibrated", action='store_true',help='Apply default calibration factors to arrays.')
parser.add_argument("--no-off", action='store_true',help='Null the off-diagonals.')
parser.add_argument("--no-prewhiten", action='store_true',help='Do not prewhiten spectra before smoothing. Use this flag for Planck.')
parser.add_argument("--overwrite", action='store_true',help='Overwrite an existing version.')
parser.add_argument("--debug", action='store_true',help='Debug plots.')
parser.add_argument("--lmax",     type=int,  default=None,help="Maximum ell.")
args = parser.parse_args()
coadd = not(args.aminusc)
nsims = args.nsims
if args.mask_patch is None: mask_patch = args.patch
else: mask_patch = args.mask_patch
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
                                      mkdir=True,overwrite=args.overwrite,mask_patch=mask_patch)
# Get data model
mask = sints.get_act_mr3_crosslinked_mask(mask_patch,
                                          version=args.mask_version,
                                          kind=args.mask_kind,
                                          season=args.season,array=args.array+"_f150",
                                          pad=args.mask_pad)

if args.debug: noise.plot(pout+"_mask",mask,grid=True)
dm = sints.models[args.model](region=mask,calibrated=args.calibrated)

# Get a NoiseGen model
if args.extract_mask is not None:
    emask = sints.get_act_mr3_crosslinked_mask(mask_patch,version=args.extract_mask,kind=args.mask_kind,season=args.season,array=args.array+"_f150")
    eshape,ewcs = emask.shape,emask.wcs
else:
    emask = mask
ngen = noise.NoiseGen(version=version,model=args.model,extract_region=emask,ncache=1,verbose=True)

# Get arrays from array

splits = dm.get_splits(season=args.season,patch=args.patch,arrays=dm.array_freqs[args.array],srcfree=True)
ivars = dm.get_splits_ivar(season=args.season,patch=args.patch,arrays=dm.array_freqs[args.array])
if args.debug: 
    noise.plot(pout+"_splits",splits)
    noise.plot(pout+"_ivars",ivars)

modlmap = splits.modlmap()
flatstring = "un" if args.do_only_filter_noise else ""
n2d_xflat = noise.get_n2d_data(splits,ivars,mask,coadd_estimator=coadd,
                               flattened=not(args.do_only_filter_noise),
                               plot_fname=pout+"_n2d_%sflat" % flatstring if args.debug else None,
                               dtype=dm.dtype)
ncomps = n2d_xflat.shape[0]
if ncomps==1: npol = 1
else: npol = 3
mask_ell = args.rlmin - args.radial_fit_annulus
del splits

radial_pairs = [(0,0),(1,1),(2,2),(3,3),(4,4),(5,5),(0,3),(3,0)] if not(args.no_prewhiten) else []
if smooth:
    n2d_xflat_smoothed = noise.smooth_ps(n2d_xflat.copy(),dfact=dfact,
                                        radial_pairs=radial_pairs,
                                        plot_fname=pout+"_n2d_%sflat_smoothed" % flatstring if args.debug else None,
                                        radial_fit_annulus = args.radial_fit_annulus,
                                        radial_fit_lmin=args.rlmin,fill_lmax=args.lmax)
else:
    n2d_xflat_smoothed = n2d_xflat.copy()
del n2d_xflat

n2d_xflat_smoothed[:,:,modlmap<mask_ell] = 0
if args.lmax is not None: n2d_xflat_smoothed[:,:,modlmap>args.lmax] = 0

if args.no_off: n2d_xflat_smoothed = noise.null_off_diagonals(n2d_xflat_smoothed)

if args.do_only_filter_noise:
    ngen.save_filter_noise(n2d_xflat_smoothed,season=args.season,patch=args.patch,array=args.array,coadd=coadd,mask_patch=mask_patch)    
    sys.exit()

covsqrt = noise.get_covsqrt(n2d_xflat_smoothed,args.covsqrt_kind)
del n2d_xflat_smoothed
ngen.save_covsqrt(covsqrt,season=args.season,patch=args.patch,array=args.array,coadd=coadd,mask_patch=mask_patch)

if nsims>0:
    bin_edges = np.arange(40,8000,40)
    p1ds = []
    for i in range(nsims):
        print("Sim %d of %d ..." % (i+1,nsims))
        with bench.show("simgen"):
            sims = ngen.generate_sim(season=args.season,patch=args.patch,array=args.array,seed=i,mask_patch=mask_patch)
            print(sims.nbytes/1024./1024./1024., " GB", sims.shape, sims.dtype)
        if args.extract_mask is not None: 
            ivars2 = enmap.extract(ivars,eshape,ewcs)
            modlmap = enmap.modlmap(eshape,ewcs)
        else:
            ivars2 = ivars

        if args.debug and i==0: noise.plot(pout+"_sims",sims)
        if not(args.no_write):
            ngen.save_sims(i,sims,args.season,args.patch,args.array,coadd=coadd,mask_patch=mask_patch)
        n2d_sim = noise.get_n2d_data(sims,ivars2,emask,coadd_estimator=coadd,
                                     flattened=False,
                                     plot_fname=pout+"_n2d_sim" if (args.debug and i==0) else None,
                                     dtype = dm.dtype)
        del sims
        cents,op1ds_sim = noise.get_p1ds(n2d_sim,modlmap,bin_edges)
        p1ds.append(op1ds_sim.copy().reshape(-1))
    p1dstats = stats.get_stats(np.array(p1ds))

    del covsqrt

    # For verification
    splits = dm.get_splits(season=args.season,patch=args.patch,arrays=dm.array_freqs[args.array],srcfree=True)

    if args.extract_mask is not None: 
        splits = enmap.extract(splits,eshape,ewcs)

    n2d_data = noise.get_n2d_data(splits,ivars2,emask,coadd_estimator=coadd,dtype=dm.dtype)
    cents,p1ds_data = noise.get_p1ds(n2d_data,modlmap,bin_edges)
    corr = noise.corrcoef(n2d_data)
    del n2d_data
    bin_edges = np.arange(40,10000,100)
    cents2,c1ds_data = noise.get_p1ds(corr,modlmap,bin_edges)
    noise.plot_corrcoeff(cents2,c1ds_data,plot_fname=pout)

    nfreqs = len(dm.array_freqs[args.array])
    noise.compare_ps(cents,p1dstats['mean'].reshape((nfreqs*npol,nfreqs*npol,cents.size)),p1ds_data,plot_fname="%s_compare" % (pout),err=p1dstats['errmean'].reshape((nfreqs*npol,nfreqs*npol,cents.size)))
