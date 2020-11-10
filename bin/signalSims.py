from __future__ import print_function
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
import time
import yaml
from enlib import bench
from pixell import aberration
import logging
from actsims import util as autil

defaults = autil.config_from_yaml("../inputParams/simsInput_v0p5_test.yaml")

import argparse
# Parse command line
parser = argparse.ArgumentParser(description='Generate lensed CMB.')
parser.add_argument("--nsims",     type=int,  default=defaults['nsims'],help="Number of sims.")
parser.add_argument("--skip-aberration", action='store_true',help='Skip aberration.')
parser.add_argument("--lmax",     type=int,  default=defaults['lmax'],help="Maxmimum multipole for lensing.")
parser.add_argument("--lmax-write",     type=int,  default=defaults['lmax_write'],help="Maximum multipole to write.")
parser.add_argument("--pix-size",     type=float,  default=defaults['pix_size'],help="Pixel width in arcminutes.")
parser.add_argument("--cmb-sets",     type=int,  nargs='+', default=defaults['cmb_sets'],help="CMB sets.")
parser.add_argument("--phi-sets",     type=int,  nargs='+', default=defaults['phi_sets'],help="phi sets.")
parser.add_argument("--input-spec",     type=str,  default=defaults['input_spec_root'],help="Input spectrum root.")
args = parser.parse_args()




def useful_info(script_filename, input_filename, parameter_dict, git_version = None, run_time = None, size = None):

    output = '\n'
    output += 'This run of %s was invoked with config file %s\n\n' % (script_filename, input_filename)

    if git_version is not None:
        output += "Git version number: %s  \n\n" % git_version
    output += "Parameters were:\n\n"
    for key in parameter_dict.keys():
        output += "%s = %s \n" % (key, str(parameter_dict[key]))
    if run_time is not None:
        output += '\nThe run took %f seconds.\n' % run_time
    if size is not None:
        output += '\nThe run was done on %i MPI processes.' % size
    return output

if __name__ == '__main__':

    try:
        with open('../inputParams/paths_local.yml') as f:
            p = yaml.safe_load(f)
    except:
        logging.error("../inputParams/paths_local.yml not found. Please copy ../inputParams/paths.yml to this file and edit with your local paths.")
        sys.exit(1)

    cmb_dir = p['data_dir']
    # logging.basicConfig(filename='/scratch/r/rbond/msyriac/example.log', level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s')


    start_time = time.clock()

    comm,rank,my_tasks = autil.distribute(args.nsims)


    shape, wcs = enmap.fullsky_geometry(args.pix_size*utils.arcmin)
    ps = powspec.read_camb_full_lens(args.input_spec + "_lenspotentialCls.dat")

    if not(args.skip_aberration):
        with bench.show("init ab"):
            ab = aberration.Aberrator(shape, wcs, modulation=None)


    #make phi totally uncorrelated with both T and E.  This is necessary due to the way that separate phi and CMB seeds were put forward in an update to the pixell library around mid-Nov 2018
    ps[0, 1:, :] = 0.
    ps[1:, 0, :] = 0.


    # if rank == 0:
    #     git_version = subprocess.check_output("git rev-parse --short HEAD", shell=True).encode().rstrip("\n")
    #     with open(cmb_dir + "/__README.txt", 'w') as f:
    #         print(useful_info(__file__, sys.argv[1], p, git_version, time.time() - start, size), file = f)


    start = time.time()

    for cmb_set in args.cmb_sets:    
        for phi_set in args.phi_sets:    
            for iii in my_tasks:
                logging.info('rank', rank, 'doing cmb_set', cmb_set, 'task' , iii, \
                    'calling lensing.rand_map', time.time() - start)


                phi_seed = seedgen.get_phi_seed(phi_set, iii)
                cmb_seed = seedgen.get_cmb_seed(cmb_set, iii)

                with bench.show("lensing"):
                    l_tqu_map, = lensing.rand_map((3,)+shape, wcs, ps,
                                                              lmax=args.lmax,
                                                              output="l",
                                                              verbose=True,
                                                              phi_seed = phi_seed,
                                                              seed = cmb_seed)



                map_list = [l_tqu_map]

                map_name_list = ['fullskyLensedUnaberratedCMB']

                if not(args.skip_aberration):


                    logging.info('doing aberration')
                    logging.info('rank', rank, 'doing cmb_set', cmb_set, 'task' , iii, \
                        'calling aberration.boost_map', time.time() - start)


                    #Note - we are aberrating and not modulating! The
                    #modulation is a frequency-dependent, so is done
                    #later.
                    with bench.show("boost"):
                        l_tqu_map_aberrated = ab.aberrate(l_tqu_map)

                    map_list += [l_tqu_map_aberrated]
                    map_name_list += ['fullskyLensedAberratedCMB']


                for mi, mmm in enumerate(map_list):
                    logging.info(iii, ' calling curvedsky.map2alm for ', map_name_list[mi])
                    alm = curvedsky.map2alm(mmm, lmax=args.lmax_write)


                    filename = cmb_dir + "/%s_alm_%s%s%05d.fits" \
                               % ( map_name_list[mi],
                                   ('cmbset%02d_' % cmb_set if 'CMB' in map_name_list[mi] else '' ) ,
                                   ('phiset%02d_' % phi_set \
                                    if (('Lensed' in map_name_list[mi]) or ('Phi' in map_name_list[mi])) \
                                    else '' ) ,
                                   iii)

                    logging.info('writing to disk, filename is', filename)

                    healpy.fitsfunc.write_alm(filename ,
                                               np.complex64(alm), overwrite = True)
