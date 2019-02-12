from __future__ import print_function
from pixell import enmap,enplot
import numpy as np
import os,sys
from actsims import datamodel, powtools
from enlib import bench
from orphics import io
import matplotlib.pyplot as plt

import argparse
# Parse command line
parser = argparse.ArgumentParser(description='Do a thing.')
parser.add_argument("season", type=str,help='Positional arg.')
parser.add_argument("array", type=str,help='Positional arg.')
parser.add_argument("patch", type=str,help='Positional arg.')
parser.add_argument("-d", "--dfact",     type=int,  default=8,help="A description.")
parser.add_argument("-c", "--coadd", action='store_true',help='Whether to use the coadd estimator.')
# parser.add_argument("-f", "--flag", action='store_true',help='A flag.')
required_args = parser.add_argument_group('Required arguments')
# required_args.add_argument("--region",     type=str,help="Region name from input/regions.yaml",required=True)
args = parser.parse_args()

dfact = (args.dfact,args.dfact)
seed = 100

bin_edges = np.arange(30,8000,40)

dm = datamodel.DataModel(args.season,args.array,args.patch)
modlmap = dm.modlmap
n2d_data = dm.get_n2d_data(dm.get_map(),coadd_estimator=args.coadd)
cents,p1ds_data = powtools.get_p1ds(n2d_data,modlmap,bin_edges)
del n2d_data
n2d_flat = dm.get_n2d_data(dm.get_map(),coadd_estimator=args.coadd,flattened=True)
n2d_flat_smoothed = powtools.smooth_ps(n2d_flat.copy(),dfact=dfact,radial_pairs=[(0,0),(1,1),(2,2),(3,3),(4,4),(5,5),(0,3)])

for corrfunc,corrlabel in zip([lambda x: x,powtools.null_off_diagonals],['allcorrs','no-off']):
    for method in ['arrayops']:
        for n2d,smlabel in zip([n2d_flat,n2d_flat_smoothed][::-1],['unsmoothed','smoothed'][::-1]):
            print(smlabel," ",method," ",corrlabel)
            covsqrt = powtools.get_covsqrt(corrfunc(n2d),method)
            sims = dm.generate_noise_sim(covsqrt,seed=seed)
            n2d_sim = dm.get_n2d_data(sims,coadd_estimator=args.coadd,flattened=False)
            cents,p1ds_sim = powtools.get_p1ds(n2d_sim,modlmap,bin_edges)
            powtools.compare_ps(cents,p1ds_sim,p1ds_data)
