from __future__ import print_function
from orphics import maps,io,cosmology,stats
from pixell import enmap
import numpy as np
import os,sys
from actsims import signal

shape,wcs = enmap.fullsky_geometry(res=np.deg2rad(4./60.))
cmb_type = 'LensedUnabberatedCMB'

model = 'planck_hybrid'
dobeam = True
apply_window = False
max_cached = 0

s1 = signal.SignalGen(cmb_type=cmb_type, dobeam=dobeam, add_foregrounds=True, apply_window=apply_window, max_cached=max_cached, model=model)

# Note freq_idx = 1 for 150 GHz
imap_143 = s1.get_signal_sim('planck', 'planck', 'planck', '143', sim_num=0,save_alm=True, save_map=False, set_idx=0,oshape=shape,owcs=wcs,fgflux="15mjy", add_poisson_srcs=False,freq_idx=1)

# Note freq_idx = 0 for 90 GHz
imap_090 = s1.get_signal_sim('planck', 'planck', 'planck', '143', sim_num=0,save_alm=True, save_map=False, set_idx=0,oshape=shape,owcs=wcs,fgflux="15mjy", add_poisson_srcs=False,freq_idx=0)


