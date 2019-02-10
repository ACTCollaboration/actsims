from orphics import maps,io,stats
from enlib import bench,array_ops
from pixell import enmap,utils
import numpy, numpy as np
import os
from tilec import covtools

seed = 100
dirpath='/home/msyriac/data/act/maps/mr3'
def load_wins(freqs,dirpath='/home/msyriac/data/act/maps/mr3',spa='s15_deep56_pa3',nsplits=4):
    """
    returns (nfreqs,nsplits,ncomp,Ny,Nx) imap
    returns (nfreqs,nsplits,1,Ny,Nx) wmap
    """
    shape,wcs = enmap.read_map_geometry("%s/%s_%s_nohwp_night_3pass_4way_set%d_map_srcfree.fits" % (dirpath,spa,freqs[0],0))
    Ny,Nx = shape[-2:]
    nfreqs = len(freqs)
    ncomp = 3
    cwmaps = np.empty((nfreqs,1,Ny,Nx))
    for ifreq,freq in enumerate(freqs):
        cwmaps[ifreq] = enmap.read_map("%s/%s_%s_nohwp_night_3pass_4way_coadd_ivar.fits" % (dirpath,spa,freq))[None]
    cwmaps = enmap.enmap(cwmaps,wcs,copy=False)
    return cwmaps

spa = 's13_deep6_pa1'
ncomp = 3
freqs = ["f150"]
region = "deep6"
cwmaps = load_wins(freqs,spa=spa)

# mtypes = [
#     # 's15_deep56_pa3_unsmoothed_pscov.fits',
#     's15_deep56_pa3_smoothed_pscov.fits']

mtypes = [
    's13_deep6_pa1_unsmoothed_pscov.fits',
    's13_deep6_pa1_smoothed_pscov.fits']

for mtype in mtypes:

    for atype in ['multipow','arrayops','multipow_noiqu','arrayops_noiqu']:
        np.random.seed(seed)
        print(mtype, atype)
        fname = mtype.replace('.fits','_%s_covsqrt.fits' % atype)
        covsqrt = enmap.read_map(fname)
        modlmap = covsqrt.modlmap()
        covsqrt[...,modlmap<2.] = 0

        covsqrt = covsqrt * np.sqrt(np.prod(covsqrt.shape[-2:]) / enmap.area(covsqrt.shape[-2:], covsqrt.wcs ))
        rmap = enmap.rand_gauss_harm((covsqrt.shape[0], covsqrt.shape[-2:][0], covsqrt.shape[-2:][1]),covsqrt.wcs)
        kmap = enmap.map_mul(covsqrt, rmap)
        spin = np.repeat([0], kmap.shape[0])
        outmaps = enmap.harm2map(kmap, iau = False, spin = spin)
        outmaps = outmaps.reshape( (len(freqs), ncomp,outmaps.shape[-2] , outmaps.shape[-1])) / np.sqrt(cwmaps)

        
        binary_percentile=10.
        bmask = maps.binary_mask(cwmaps,threshold = np.percentile(cwmaps,binary_percentile))
        for i in range(len(freqs)):
            outmaps[i,:,bmask[i,0]==0] = 0
        # io.plot_img(outmaps,lim=300)
        freq = freqs[0]
        imap = enmap.read_map("%s/%s_%s_nohwp_night_3pass_4way_coadd_map_srcfree.fits" % (dirpath,spa,freq))
        io.plot_img(imap,lim=300)

        enmap.write_map(fname.replace(".fits","_sim_seed_%d.fits"%seed),enmap.enmap(outmaps,rmap.wcs))

