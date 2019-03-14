from actsims import signal, noise
from pixell import enmap,wcsutils
import numpy as np


"""

cmbSeedInd = 0
fgSeedInd = 1
phiSeedInd = 2
noiseSeedInd = 3

phiSeed = (0, 0, 2, i)
cmbSeed = (set, 0, 0, i)
fgseed = (set, 0, 1, i)

"""


class SimGen(object):
    def __init__(self, version, model="act_mr3", cmb_type='LensedCMB', dobeam=True, add_foregrounds=True, apply_window=True, max_cached=1):
        
        """
        version: The version identifier for the filename of covsqrts on disk
        model: The name of an implemented soapack datamodel
        extract_region: An optional map whose footprint on to which the sims are made
        extract_region_shape: Instead of passing a map for extract_region, one can pass its shape and wcs
        extract_region_wcs: Instead of passing a map for extract_region, one can pass its shape and wcs
        max_cached: The maximum number of cached sim/or alms
        """
 
        
        self.noise_gen  = noise.NoiseGen(version=version,model=model,ncache=max_cached,verbose=False)
        self.signal_gen = signal.SignalGen(cmb_type=cmb_type, dobeam=dobeam, add_foregrounds=add_foregrounds, apply_window=apply_window, max_cached=max_cached, model=model)


    ## Note: It will be probably better to use a decorator to delegate. For now, it does it explicetly

    def get_signal(self, season, patch, array, freq, sim_num, save_alm=True, save_map=False, set_idx=0,oshape=None,owcs=None):
        # return cmb+fg sim
        return self.signal_gen.get_signal_sim(season, patch, array, freq, set_idx, sim_num, save_alm, save_map,oshape=oshape,owcs=owcs)

    def get_cmb(self, season, array, patch, freq, sim_num, save_alm=False, set_idx=0,oshape=None,owcs=None):
        return self.signal_gen.get_cmb_sim(season, array, patch, freq, set_idx, sim_num, save_alm,oshape=oshape,owcs=owcs)

    def get_fg(self, season, array, patch, freq, sim_num, save_alm=False, set_idx=0,oshape=None,owcs=None):
        return self.signal_gen.get_fg_sim(season, array, patch, freq, set_idx, sim_num, save_alm,oshape=oshape,owcs=owcs)

    def get_noise(self, season=None,patch=None,array=None, sim_num=None,mask_patch=None,set_idx=0,apply_ivar=True):
        # indexing is slighly different for signal and noise sim code ..
        # array = None if (array is not None) and (freq is not None) else '{}_{}'.format(array, freq)   
        seed = (set_idx, 0, 3, sim_num) + self.noise_gen.dm.get_noise_sim_seed(season,patch,array)
        print(seed)
        return self.noise_gen.generate_sim(season,patch,array,seed=seed,mask_patch=mask_patch,apply_ivar=apply_ivar)

    def get_sim(self,season,patch,array,sim_num, save_alm=True, save_map=False, set_idx=0,mask_patch=None):
        shape,wcs = self.noise_gen.load_covsqrt(season,patch,array,coadd=True,mask_patch=mask_patch,get_geometry=True)
        # (nfreqs,nsplits,npol,Ny,Nx)
        noises,ivars = self.get_noise(season=season,patch=patch,array=array, sim_num=sim_num,mask_patch=mask_patch,set_idx=set_idx,apply_ivar=False)
        from orphics import io
        io.hplot(ivars,"ivars")
        afreqs = self.noise_gen.dm.array_freqs[array]
        signals = []
        for afreq in afreqs:
            imap = self.get_signal(season, patch, array, afreq.split('_')[1], sim_num, save_alm=save_alm, save_map=save_map, set_idx=set_idx,oshape=shape,owcs=wcs)
            print(afreq,imap.shape,imap.wcs)
            signals.append(imap)
        owcs = imap.wcs
        # (nfreqs,npol,Ny,Nx)
        signals = enmap.enmap(np.stack(signals),owcs)[:,None,...]
        del imap
        assert wcsutils.equal(wcs,noises.wcs)
        assert wcsutils.equal(owcs,noises.wcs)
        return noise.apply_ivar_window(signals + noises,ivars)
