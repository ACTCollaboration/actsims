"""
This script was used to discover the usefulness of the difference estimator and of the noise PS smoothing.
"""

from __future__ import print_function
from pixell import enmap,enplot
import numpy as np
import os,sys
from actsims import datamodel, powtools,utils
from enlib import bench
from orphics import io,stats
import matplotlib.pyplot as plt
from tilec import covtools
import argparse

# Parse command line
parser = argparse.ArgumentParser(description='Make covsqrt, generate some test sims, make verification plots.')
parser.add_argument("season", type=str,help='Season')
parser.add_argument("array", type=str,help='Array')
parser.add_argument("patch", type=str,help='Patch')
parser.add_argument("-n", "--nsims",     type=int,  default=10,help="Number of sims.")
parser.add_argument("-d", "--dfact",     type=int,  default=8,help="Downsample factor.")
parser.add_argument("-a", "--aminusc", action='store_true',help='Whether to use the auto minus cross estimator.')
parser.add_argument("--no-write", action='store_true',help='Do not write any FITS to disk.')
args = parser.parse_args()
coadd = not(args.aminusc)
dfact = (args.dfact,args.dfact)
nsims = args.nsims

bin_edges = np.arange(30,8000,40)

dm = datamodel.DataModel(args.season,args.array,args.patch)
pout = "%s%s_%s_%s_coadd_est_%s" % (datamodel.pout ,args.season,args.array,args.patch,coadd)
sout = "%s%s_%s_%s_coadd_est_%s" % (datamodel.paths['save'] ,args.season,args.array,args.patch,coadd)

modlmap = dm.modlmap

n2d_data = dm.get_n2d_data(dm.get_map(),coadd_estimator=coadd)#,plot_fname=pout+"_n2d_data")
cents,p1ds_data = powtools.get_p1ds(n2d_data,modlmap,bin_edges)
del n2d_data

n2d_flat = dm.get_n2d_data(dm.get_map(),coadd_estimator=coadd,flattened=True)#,plot_fname=pout+"_n2d_data_flat")
n2d_flat_smoothed = powtools.smooth_ps(n2d_flat.copy(),dfact=dfact,radial_pairs=[(0,0),(1,1),(2,2),(3,3),(4,4),(5,5),(0,3),(3,0)])#,plot_fname=pout+"_n2d_data_flat_smoothed")

for corrfunc,corrlabel in zip([lambda x: x,powtools.null_off_diagonals],['allcorrs','no-off']):
    for method in ['arrayops']:
        for n2d,smlabel in zip([n2d_flat,n2d_flat_smoothed][::-1],['unsmoothed','smoothed'][::-1]):
            print(smlabel," ",method," ",corrlabel)
            covsqrt = powtools.get_covsqrt(corrfunc(n2d),method)
            if smlabel=='smoothed' and corrlabel=='allcorrs' and method=='arrayops' and not(args.no_write):
                enmap.write_map("%s_covsqrt_%s_%s_%s.fits" % (sout,corrlabel,method,smlabel) ,covsqrt)
            p1ds = []
            for i in range(nsims):
                print("Sim %d of %d ..." % (i+1,nsims))
                with bench.show("simgen"):
                    sims = dm.generate_noise_sim(covsqrt,seed=i)
                if smlabel=='smoothed' and corrlabel=='allcorrs' and method=='arrayops' and not(args.no_write):
                    enmap.write_map("%s_trial_sim_%s_%s_%s_seed_%d.fits" % (sout,corrlabel,method,smlabel,i) ,sims)
                n2d_sim = dm.get_n2d_data(sims,coadd_estimator=coadd,flattened=False)#,plot_fname=pout+"_n2d_sims")
                del sims
                cents,op1ds_sim = powtools.get_p1ds(n2d_sim,modlmap,bin_edges)
                p1ds.append(op1ds_sim.copy().reshape(-1))
            p1dstats = stats.get_stats(np.array(p1ds))
            powtools.compare_ps(cents,p1dstats['mean'].reshape((dm.nfreqs*3,dm.nfreqs*3,cents.size)),p1ds_data,plot_fname="%s%s_%s_%s_compare" % (pout,corrlabel,method,smlabel),err=p1dstats['errmean'].reshape((dm.nfreqs*3,dm.nfreqs*3,cents.size)))
