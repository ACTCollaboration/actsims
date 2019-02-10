import pitas
from actsims import simTools
from orphics import maps, io
from pixell import enmap, curvedsky
import sys, os
import numpy as np
from itertools import product
from actsims import simTools
import pickle

region = "deep6"

xlink = enmap.read_map("/home/msyriac/data/act/maps/steve/%s_mask_run_180323_master_apo_w0.fits" % region)

output_dir     = '/global/homes/d/dwhan89/shared/outbox/cori/121018_coadded_realstic'
output_path    = lambda x: os.path.join(output_dir, x)

patchlist      = ['s13_deep6_pa1_f150']

mask           = xlink
modlmap        = mask.modlmap()

lmax           = 4000
bin_edges      = pitas.util.get_default_bin_edges(lmax)
ret_dl         = True

PITAS    = pitas.power.PITAS("coadded_noise_test%d"%lmax, mask, mask, bin_edges, lmax, transfer=None, overwrite=False)

def add_cmb_spec(st, idx, tmap1, emap1, bmap1, tmap2=None, emap2=None, bmap2=None, key_postfix='', overwrite=True, ret_dl=ret_dl):
    key_temp  = 'dl{}_%s' %key_postfix
    global PITAS
    lbin = PITAS.bin_center
    dltt, dlte, dlee, dlbb = (None, None, None, None)
    if not st.has_data(key_temp.format('bb'), idx) or overwrite:
        if tmap2 is None:
            lbin, dltt = PITAS.get_power(tmap1, tmap1, polcomb='TT', ret_dl=ret_dl)
            lbin, dlte = PITAS.get_power(tmap1, emap1, polcomb='TE', ret_dl=ret_dl)
            lbin, dlee, dleb, dlbb = PITAS.get_power_pureeb(emap1, bmap1, ret_dl=ret_dl)
        else:
            lbin, dltt = PITAS.get_power(tmap1, tmap2, polcomb='TT', ret_dl=ret_dl)
            lbin, dlte1 = PITAS.get_power(tmap1, emap2, polcomb='TE', ret_dl=ret_dl)
            lbin, dlte2 = PITAS.get_power(tmap2, emap1, polcomb='TE', ret_dl=ret_dl)
            dlte = (dlte1+dlte2)/2.
            lbin, dlee, dleb, dlbb = PITAS.get_power_pureeb(emap1, bmap1, emap2, bmap2, ret_dl=ret_dl)

        st.add_data(key_temp.format('tt'), idx, dltt)
        st.add_data(key_temp.format('te'), idx, dlte)
        st.add_data(key_temp.format('ee'), idx, dlee)
        st.add_data(key_temp.format('bb'), idx, dlbb)
    else:
        _, dltt, dlte, dlee, dlbb = get_cmb_spec(st, idx, key_postfix)

    return (lbin, dltt, dlte, dlee, dlbb)

def get_noise_sims(fname, idx,calF=1., calFPol=1.):
    seed = 100
    noise = enmap.read_map(fname.replace(".fits","_sim_seed_%d.fits"%seed))

    noise[:, 0, :, :] = noise[:, 0, :, :] * np.sqrt(calF)
    noise[:, 1, :, :] = noise[:, 1, :, :] * np.sqrt(calFPol)
    noise[:, 2, :, :] = noise[:, 2, :, :] * np.sqrt(calFPol)

    ret = {}
    ret['i'] = noise[idx, 0, :, :].copy()
    ret['q'] = noise[idx, 1, :, :].copy()
    ret['u'] = noise[idx, 2, :, :].copy()
    return ret


def tqu2teb(tmap, qmap, umap, lmax=8000):
    emap = bmap = None

    tqu     = np.zeros((3,)+tmap.shape)
    tqu[0], tqu[1], tqu[2] = (tmap, qmap, umap)
    tqu     = enmap.enmap(tqu, tmap.wcs)
    alm     = curvedsky.map2alm(tqu, lmax=lmax)

    teb     = curvedsky.alm2map(alm[:,None], tqu.copy()[:,None], spin=0)[:,0]
    del tqu

    return (tmap, teb[1], teb[2])

mtypes = [
    's13_deep6_pa1_unsmoothed_pscov.fits',
    's13_deep6_pa1_smoothed_pscov.fits']
freqs = ["f150"]

from orphics import io
pl = io.Plotter()
for mtype in mtypes:
    for atype in ['multipow','arrayops','multipow_noiqu','arrayops_noiqu']:
        for idx in range(len(freqs)):
            
            fname = mtype.replace('.fits','_%s_covsqrt.fits' % atype)
            temp = get_noise_sims(fname,idx)
            for cmb_idx in ['i', 'q', 'u']:
                temp[cmb_idx]  = temp[cmb_idx]*mask

            _, temp['e'], temp['b'] = tqu2teb(temp['i'].copy(), temp['q'], temp['u'])
            tmap1 = temp['i']
            emap1 = temp['e']
            bmap1 = temp['b']
            
            lbin, dltt = PITAS.get_power(tmap1, tmap1, polcomb='TT', ret_dl=ret_dl)
            # lbin, dlte = PITAS.get_power(tmap1, emap1, polcomb='TE', ret_dl=ret_dl)
            # lbin, dlee, dleb, dlbb = PITAS.get_power_pureeb(emap1, bmap1, ret_dl=ret_dl)

            pl.add(lbin,dltt,label=mtype+atype)
pl.done()

