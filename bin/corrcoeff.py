"""
This script reads saved simulations and re-makes verification plots.
"""

from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
from pixell import enmap,enplot
import numpy as np
import os,sys
from actsims import noise as datamodel
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

bin_edges = np.arange(30,10000,100)

dm = datamodel.NoiseModel(args.season,args.array,args.patch)
pout = "%s%s_%s_%s_coadd_est_%s" % (datamodel.pout ,args.season,args.array,args.patch,True)
sout = "%s%s_%s_%s_coadd_est_%s" % (datamodel.paths['save'] ,args.season,args.array,args.patch,True)

modlmap = dm.modlmap

n2d_data = dm.get_n2d_data(dm.get_map(),coadd_estimator=True)
corr = datamodel.corrcoef(n2d_data)
cents,c1ds_data = datamodel.get_p1ds(corr,modlmap,bin_edges)

dpi = 300
ncomps = c1ds_data.shape[0]
if ncomps==3:
    pols = ['150-I','150-Q','150-U']
elif ncomps==6:
    pols = ['90-I','90-Q','90-U','150-I','150-Q','150-U']

pl = io.Plotter(xlabel = "$\\ell$", ylabel = "$N_{XY}/\\sqrt{N_{XX}N_{YY}}$",xyscale='linlin')
for i in range(c1ds_data.shape[0]):
    for j in range(i+1,c1ds_data.shape[0]):
        polstring = "%s x %s" % (pols[i],pols[j])
        pl.add(cents,c1ds_data[i,j],label=polstring)
pl._ax.set_xlim(30,10000)
pl._ax.set_ylim(-0.3,0.3)
pl.hline(y=0.05)
pl.hline(y=-0.05)
pl.hline(y=0.01,ls='-.')
pl.hline(y=-0.01,ls='-.')
pl.legend(loc='center left', bbox_to_anchor=(1, 0.5))
pl.hline(y=0,ls='-')
pl.vline(x=500)
pl.done("%s_corrcoeff.png" % (pout), dpi=dpi)


