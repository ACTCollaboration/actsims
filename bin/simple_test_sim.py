from __future__ import print_function
from orphics import io,stats
from pixell import enmap
import numpy as np
import os,sys


def rednoise(ells,rms_noise,lknee=0.,alpha=1.):
    """Atmospheric noise model
    rms_noise in muK-arcmin
    [(lknee/ells)^(-alpha) + 1] * rms_noise**2
    """
    atm_factor = (lknee*np.nan_to_num(1./ells))**(-alpha) if lknee>1.e-3 else 0.
    rms = rms_noise * (1./60.)*(np.pi/180.)
    wnoise = ells*0.+rms**2.
    return (atm_factor+1.)*wnoise

def get_covsqrt(ps,method="arrayops"):
    if method=="multipow":
        covsq = enmap.multi_pow(ps.copy(),0.5)
    elif method=="arrayops":
        from enlib import array_ops
        covsq = array_ops.eigpow(ps.copy(),0.5,axes=[0,1])
    covsq[:,:,ps.modlmap()<2] = 0
    assert np.all(np.isfinite(covsq))
    return covsq

patch = "deep6"


masks = {}
masks['new'] = enmap.read_map("%s_mr3c_20181012_190203_master_apo_w0.fits" % patch)
masks['old'] = enmap.read_map("../180323/%s_mask_run_180323_master_apo_w0.fits" % patch)

pl = io.Plotter(xyscale='linlin')
for method in ['arrayops','multipow']:
    for amask in ['old','new']:

        mask = masks[amask]
        shape,wcs = mask.shape,mask.wcs
        modlmap = enmap.modlmap(shape,wcs)
        Ny,Nx = shape[-2:]

        n2d = rednoise(modlmap,rms_noise=20.,lknee=3000.,alpha=-4.5)
        n2d[modlmap<100] = 0

        bin_edges = np.arange(100,8000,40)
        binner = stats.bin2D(modlmap,bin_edges)
        cents,n1d = binner.bin(n2d)

        covsqrt = get_covsqrt(n2d[None,None].copy())
        rmap = enmap.rand_gauss_harm((1, Ny, Nx),covsqrt.wcs) 
        kmap = enmap.map_mul(covsqrt, rmap)
        imap = enmap.ifft(kmap, normalize="phys").real

        kmap = enmap.fft(imap*mask,normalize="phys")
        p2d = np.real(kmap*np.conjugate(kmap))/np.mean(mask**2.)


        cents,p1d = binner.bin(p2d)

        pl.add(cents,p1d/n1d,label=method+amask)
pl.hline(y=1)
pl.done()
        
# pl = io.Plotter(xyscale='linlog')
# pl.add(cents,p1d)
# pl.add(cents,n1d)
# pl.hline(y=1)
# pl.done()

