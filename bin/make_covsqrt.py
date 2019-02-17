"""
This script can be used to make a covsqrt and a few trial sims.
"""
from __future__ import print_function
from pixell import enmap,enplot
import numpy as np
import os,sys
from actsims import noise
from soapack import covsqrt
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
parser.add_argument("--debug", action='store_true',help='Debug plots.')
args = parser.parse_args()
coadd = not(args.aminusc)
dfact = (args.dfact,args.dfact)
nsims = args.nsims

bin_edges = np.arange(30,8000,40)

dm = datamodel.DataModel(args.season,args.array,args.patch)
pout = "%s%s_%s_%s_coadd_est_%s" % (datamodel.pout ,args.season,args.array,args.patch,coadd)
sout = "%s%s_%s_%s_coadd_est_%s" % (datamodel.paths['save'] ,args.season,args.array,args.patch,coadd)

modlmap = dm.modlmap

maps = dm.get_map()
n2d_flat = dm.get_n2d_data(maps,coadd_estimator=coadd,flattened=True,plot_fname=pout+"_n2d_flat" if args.debug else None)
del maps
n2d_flat_smoothed = powtools.smooth_ps(n2d_flat.copy(),dfact=dfact,radial_pairs=[(0,0),(1,1),(2,2),(3,3),(4,4),(5,5),(0,3),(3,0)],plot_fname=pout+"_n2d_flat_smoothed" if args.debug else None)
del n2d_flat

covsqrt = powtools.get_covsqrt(n2d_flat_smoothed,"arrayops")
enmap.write_map("%s_covsqrt.fits" % (sout) ,covsqrt)

p1ds = []
for i in range(nsims):
    print("Sim %d of %d ..." % (i+1,nsims))
    with bench.show("simgen"):
        sims = dm.generate_noise_sim(covsqrt,seed=i)#,binary_percentile=50. if (args.array=='pa1' and args.patch=='deep8') else 10.)
    enmap.write_map("%s_trial_sim_seed_%d.fits" % (sout,i) ,sims)
    n2d_sim = dm.get_n2d_data(sims,coadd_estimator=coadd,flattened=False,plot_fname=pout+"_n2d_sim" if args.debug else None)
    del sims
    cents,op1ds_sim = powtools.get_p1ds(n2d_sim,modlmap,bin_edges)
    p1ds.append(op1ds_sim.copy().reshape(-1))
p1dstats = stats.get_stats(np.array(p1ds))

del covsqrt

# For verification
n2d_data = dm.get_n2d_data(dm.get_map(),coadd_estimator=coadd)
cents,p1ds_data = powtools.get_p1ds(n2d_data,modlmap,bin_edges)
del n2d_data

powtools.compare_ps(cents,p1dstats['mean'].reshape((dm.nfreqs*3,dm.nfreqs*3,cents.size)),p1ds_data,plot_fname="%s_compare" % (pout),err=p1dstats['errmean'].reshape((dm.nfreqs*3,dm.nfreqs*3,cents.size)))
