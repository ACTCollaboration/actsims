"""
This script reads saved simulations and re-makes verification plots.
"""

from __future__ import print_function
from pixell import enmap,enplot
import numpy as np
import os,sys
from actsims import noise
from soapack import 
from enlib import bench
from orphics import io,stats
import matplotlib.pyplot as plt
from tilec import covtools
import argparse
import glob

# Parse command line
parser = argparse.ArgumentParser(description='Make verification plots.')
parser.add_argument("season", type=str,help='Season')
parser.add_argument("array", type=str,help='Array')
parser.add_argument("patch", type=str,help='Patch')
args = parser.parse_args()

bin_edges = np.arange(30,10000,20)

dm = datamodel.NoiseModel(args.season,args.array,args.patch)
pout = "%s%s_%s_%s_coadd_est_%s" % (datamodel.pout ,args.season,args.array,args.patch,True)
sout = "%s%s_%s_%s_coadd_est_%s" % (datamodel.paths['save'] ,args.season,args.array,args.patch,True)

modlmap = dm.modlmap

nsims = 2 #len(glob.glob("%s_trial_sim_seed_*.fits" % sout))

p1ds = []
c1ds = []
for i in range(nsims):
    sims = enmap.read_map("%s_trial_sim_seed_%d.fits" % (sout,i))
    n2d_sim = dm.get_n2d_data(sims,coadd_estimator=True,flattened=False)
    del sims
    cents,op1ds_sim = datamodel.get_p1ds(n2d_sim,modlmap,bin_edges)
    corr = datamodel.corrcoef(n2d_sim)
    del n2d_sim
    cents,oc1ds_sim = datamodel.get_p1ds(corr,modlmap,bin_edges)
    del corr
    p1ds.append(op1ds_sim.copy().reshape(-1))
    c1ds.append(oc1ds_sim.copy().reshape(-1))
    print("Done with ", i+1, " / ", nsims , " ...")
p1dstats = stats.get_stats(np.array(p1ds))
c1dstats = stats.get_stats(np.array(c1ds))

n2d_data = dm.get_n2d_data(dm.get_map(),coadd_estimator=True)
cents,p1ds_data = datamodel.get_p1ds(n2d_data,modlmap,bin_edges)
del n2d_data

datamodel.compare_ps(cents,p1dstats['mean'].reshape((dm.nfreqs*3,dm.nfreqs*3,cents.size)),p1ds_data,plot_fname="%s_better_compare" % (pout),err=p1dstats['errmean'].reshape((dm.nfreqs*3,dm.nfreqs*3,cents.size)))

pl = io.Plotter(xlabel = "$\\ell$", ylabel = "$C$",xyscale='linlin')

