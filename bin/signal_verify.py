from __future__ import print_function
from orphics import maps,io,cosmology,stats
from pixell import enmap
import numpy as np
import os,sys
import healpy as hp
import logging
from actsims import util as autil
import yaml
from plottools.plotutils import colorscale, update_rcParams
update_rcParams()

if __name__ == '__main__':

    try:
        with open('../inputParams/paths_local.yml') as f:
            p = yaml.safe_load(f)
    except:
        logging.error("../inputParams/paths_local.yml not found. Please copy ../inputParams/paths.yml to this file and edit with your local paths.")
        sys.exit(1)

    cmb_dir = p['data_dir']

    fells,ftt = np.loadtxt('../data/fg.dat',usecols=[0,1],unpack=True)
    ftt = ftt/fells/(fells+1.)*2.*np.pi


    fnames = [f'{cmb_dir}/fullskyLensedUnabberatedCMB_alm_cmbset00_phiset00_00000.fits',f'{cmb_dir}/../../alex/v0.4/fullskyLensedUnabberatedCMB_alm_set00_00000.fits']
    versions = ['v0.5','v0.4']

    pl = io.Plotter('rCell',ylabel='$(C_L-C^{\\rm theory}_L) / (C^{\\rm theory}_L + C^{\\rm fg}_L)$')
    pl2 = io.Plotter('Dell')

    theory = cosmology.default_theory()
    ls = np.arange(10000)
    ftt = maps.interp(fells,ftt)(ls)
    cltt = theory.lCl('TT',ls)

    for i,fname in enumerate(fnames):
        alm = hp.read_alm(fname).astype(np.complex128)
        cls = hp.alm2cl(alm)
        cls = maps.interp(np.arange(cls.size),cls)(ls)

        pl._ax.scatter(ls,(cls-cltt)/(cltt+ftt),s=1)
        pl.add([],[],label=versions[i])
        pl2.add(ls,cls,label=versions[i])

    pl._ax.set_ylim(-0.02,0.02)
    pl.hline(y=0)
    pl._ax.set_xlim(0,10000)
    pl.legend(loc='upper right')
    pl.vline(x=8250)
    pl.done('rcl.png')

    pl2.add(ls,ftt,ls=':',label='Foregrounds')
    pl2.add(ls,cltt,label='Lensed theory')
    pl2._ax.set_xlim(0,10000)
    pl2.vline(x=8250)
    pl2.done('ccl.png')
