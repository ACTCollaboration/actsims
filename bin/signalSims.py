from __future__ import print_function
import logging
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from actsims.util import seed_tracker as seedgen
from pixell import enmap, utils , lensing
from pixell import powspec, curvedsky
import numpy as np
import healpy
import subprocess, sys
from mpi4py import MPI
import yaml
from enlib import bench
from pixell import aberration
from actsims import util as autil
import os,sys,json,shutil

defaults = autil.config_from_yaml("../inputParams/defaults_lcmb.yaml")
try:
    with open('../inputParams/paths_local.yml') as f:
        p = yaml.safe_load(f)
except:
    logging.error("../inputParams/paths_local.yml not found. Please copy ../inputParams/paths.yml to this file and edit with your local paths.")
    sys.exit(1)


import argparse
# Parse command line
parser = argparse.ArgumentParser(description='Generate lensed CMB.',
                                 epilog='Example: python signalSims.py --skip-aberration --cmb-phi-sets [[0,0],[1,0]] --nsims [2000,500]')
parser.add_argument("--nsims",     type=json.loads,  default=defaults['nsims'],help="List of number of sims for each cmb_phi_set.")
parser.add_argument("--start-index",     type=json.loads,  default=defaults['start_index'],help="List of start indices for each cmb_phi_set.")
parser.add_argument("--skip-aberration", action='store_true',help='Skip aberration.')
parser.add_argument("--only-show-njobs", action='store_true',help='Do not simulate; just show total number of jobs.')
parser.add_argument("--lmax",     type=int,  default=defaults['lmax'],help="Maxmimum multipole for lensing.")
parser.add_argument("--lmax-write",     type=int,  default=defaults['lmax_write'],help="Maximum multipole to write.")
parser.add_argument("--pix-size",     type=float,  default=defaults['pix_size'],help="Pixel width in arcminutes.")
parser.add_argument("--cmb-phi-sets",     type=json.loads, default=defaults['cmb_phi_sets'],help="List of cmb-phi-set index pairs. e.g. [[0,0],[1,0]]")
parser.add_argument("--input-spec",     type=str,  default=defaults['input_spec_root'],help="Input spectrum root.")
parser.add_argument("--output-dir",     type=str,  default=p['data_dir'],help="Output directory.")
args = parser.parse_args()

cmb_dir = args.output_dir

# Save args for later reference
with open(f'{cmb_dir}/args.yml','w') as f:
    f.write(yaml.dump(vars(args)))

# Create a list of jobs and MPI distribute
jobs = []
assert len(args.nsims)==len(args.cmb_phi_sets)==len(args.start_index)
for cpset,nsim,sindex in zip(args.cmb_phi_sets,args.nsims,args.start_index):
    cmb_set,phi_set = cpset
    for i in range(nsim):
        jobs.append( (cmb_set,phi_set,i+sindex) )

nsims = len(jobs)
if args.only_show_njobs:
    print(f"Number of jobs: {nsims}")
    sys.exit()
comm,rank,my_tasks = autil.distribute(nsims)

# Start a logger
if rank==0:
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(filename=f'{cmb_dir}/log.txt', level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s',filemode='w')

# Initialize geometry
shape, wcs = enmap.fullsky_geometry(args.pix_size*utils.arcmin)

# Load theory file and save for later reference
ps = powspec.read_camb_full_lens(args.input_spec + "_lenspotentialCls.dat")
shutil.copyfile(args.input_spec + "_lenspotentialCls.dat",f'{cmb_dir}/lenspotentialCls.dat')
#make phi totally uncorrelated with both T and E.  This is necessary due to the way that separate phi and CMB seeds were put forward in an update to the pixell library around mid-Nov 2018
ps[0, 1:, :] = 0.
ps[1:, 0, :] = 0.

# Initialize aberrator
if not(args.skip_aberration):
    with bench.mark("init ab"):
        ab = aberration.Aberrator(shape, wcs, modulation=None)
    if rank==0:
        logging.info(f'BENCH:\n{bench.stats}')

# Log package info
if rank == 0:
    logging.info("Saving package info...")
    logging.info(autil.pretty_info(autil.get_info(path=os.path.realpath(__file__))))
    logging.info(autil.pretty_info(autil.get_info(package='pixell')))


# Loop over tasks
for j,task in enumerate(my_tasks):

    # Get CMB and Phi seeds
    cmb_set,phi_set,iii = jobs[task]
    cmb_seed = seedgen.get_cmb_seed(cmb_set, iii) 
    phi_seed = seedgen.get_phi_seed(phi_set, iii)
    logging.info(f'rank {rank}, task {task}, doing cmb_set {cmb_set}, phi_set {phi_set}, iteration {iii}')

    # Make lensed map
    with bench.mark("lensing"):
        l_tqu_map, = lensing.rand_map((3,)+shape, wcs, ps,
                                      lmax=args.lmax,
                                      output="l",
                                      phi_seed = phi_seed,
                                      seed = cmb_seed,
                                      verbose = (True if rank==0 else False))
    if rank==0: logging.info(f'BENCH:\n{bench.stats}')

    map_list = [l_tqu_map]
    map_name_list = ['fullskyLensedUnaberratedCMB']

    if not(args.skip_aberration):

        if rank==0:
            logging.info('doing aberration')

        #Note - we are aberrating and not modulating! The
        #modulation is a frequency-dependent, so is done
        #later.
        with bench.mark("boost"):
            l_tqu_map_aberrated = ab.aberrate(l_tqu_map)
        if rank==0:
            logging.info(f'BENCH:\n{bench.stats}')


        map_list += [l_tqu_map_aberrated]
        map_name_list += ['fullskyLensedAberratedCMB']


    for mi, mmm in enumerate(map_list):
        if rank==0:
            logging.info(f'curvedsky.map2alm for {map_name_list[mi]}')
        alm = curvedsky.map2alm(mmm, lmax=args.lmax_write)
        filename = cmb_dir + f"/{map_name_list[mi]}_alm_cmb_set_{cmb_set:02d}_phi_set_{phi_set:02d}_{iii:05d}.fits"

        if rank==0:
            logging.info(f'writing to disk, filename is {filename}')

        healpy.fitsfunc.write_alm(filename ,
                                   np.complex64(alm), overwrite = True)

    if rank==0:
        logging.info(f'rank 0: {(j+1.)*100./len(my_tasks)} % done...')
