import numpy as np
import os,sys
from pixell import enmap,enplot
import powtools
from enlib import bench

map_root = "/home/msyriac/data/act/maps/mr3/"
mask_root = "/home/msyriac/data/act/maps/steve/"

plot = lambda imap,dg=8 : enplot.show(enplot.plot(enmap.downgrade(imap,dg),grid=False))

class DataModel(object):
    def __init__(self,season,array,patch):
        self.mask_a = enmap.read_map("%s%s_mask_run_180323_master_apo_w0.fits" % (mask_root,patch))
        self.shape,self.wcs = self.mask_a.shape,self.mask_a.wcs
        self.modlmap = enmap.modlmap(self.shape,self.wcs)
        self.region = patch
        self.array = array
        self.season = season
        self.freqs = {'pa1':['f150'],'pa2':['f150'],'pa3':['f090','f150']}[array]
        self.nfreqs = len(self.freqs)
        self.power = powtools.Power(self.shape,self.wcs,mc=False)
        
        
    def get_inv_var(self):
        orets = []
        for freq in self.freqs:
            rets = []
            for k in range(4):
                pref = '_'.join([self.season,self.region,self.array,freq])
                rets.append( enmap.read_map("%s%s_nohwp_night_3pass_4way_set%d_ivar.fits" % (map_root,pref,k))[None] )
                # DEBUGGING rets.append( enmap.read_map("%s%s_nohwp_night_3pass_4way_set0_ivar.fits" % (map_root,pref))[None]*0+1. )

            orets.append(np.stack(rets))
        return enmap.enmap(np.stack(orets),self.wcs)
    
    def get_map(self):
        orets = []
        for freq in self.freqs:
            rets = []
            for k in range(4):
                pref = '_'.join([self.season,self.region,self.array,freq])
                rets.append( enmap.read_map("%s%s_nohwp_night_3pass_4way_set%d_map_srcfree.fits" % (map_root,pref,k)) )
            orets.append(np.stack(rets))
        return enmap.enmap(np.stack(orets),self.wcs)


    def get_n2d_data(self,splits,coadd_estimator=False,flattened=False,show=False):
        ivars = self.get_inv_var()
        if coadd_estimator:
            coadd,_ = powtools.get_coadd(splits,ivars,axis=1)
            data  = splits - coadd[:,None,...]
            del coadd
        else:
            data = splits
        if flattened:
            ffts = enmap.fft(data*self.mask_a*np.sqrt(ivars),normalize="phys")
            wmaps = self.mask_a + enmap.zeros(ffts.shape)
            del ivars, data, splits
        else:
            ffts = enmap.fft(data*self.mask_a*ivars,normalize="phys")
            wmaps = ivars * self.mask_a
            del ivars, data, splits
        return self.power.get_n2d(ffts,wmaps,show=show,coadd_estimator=coadd_estimator)


    def generate_noise_sim(self,icovsqrt,binary_percentile=10.,seed=None):
        if isinstance(seed,int): seed = [seed]

        modlmap = enmap.modlmap(self.shape,self.wcs)
        Ny,Nx = self.shape[-2:]
        nfreqs = self.nfreqs
        ncomps = nfreqs * 3
        wmaps = self.get_inv_var()
        wcs = wmaps.wcs

        nsplits = wmaps.shape[1]

        # Old way with loop
        covsqrt = icovsqrt 
        kmap = []
        for i in range(nsplits):
            if seed is None:
                np.random.seed(None)
            else:
                np.random.seed(seed+[i])
            rmap = enmap.rand_gauss_harm((ncomps, Ny, Nx),covsqrt.wcs) 
            kmap.append( enmap.map_mul(covsqrt, rmap) )
        kmap = enmap.enmap(np.stack(kmap),self.wcs)
        outmaps = enmap.ifft(kmap, normalize="phys").real

        # Need to test this more
        # covsqrt = icovsqrt 
        # kmap = []
        # np.random.seed(seed)
        # rmap = enmap.rand_gauss_harm((nsplits,ncomps,Ny, Nx),covsqrt.wcs)
        # kmap = enmap.samewcs(np.einsum("abyx,cbyx->cayx", covsqrt, rmap),rmap)
        # outmaps = enmap.ifft(kmap, normalize="phys").real

        
        fmaps = []
        for ifreq in range(nfreqs):
            omaps = outmaps[:,ifreq*3:(ifreq+1)*3,...].copy() / np.sqrt(wmaps[ifreq,...]) *np.sqrt(nsplits)
            fmaps.append(omaps.copy())
        fmaps = enmap.enmap(np.stack(fmaps),self.wcs)
        del omaps,outmaps

        # Sanitize
        for ifreq in range(nfreqs):
            for isplit in range(nsplits):
                win = wmaps[ifreq,isplit,0,...]
                bmask = powtools.binary_mask(win,threshold = np.percentile(win,binary_percentile))
                fmaps[...,bmask==0] = 0

        return fmaps

    
