from actsims import signal, noise
from pixell import enmap,wcsutils
import numpy as np


"""

cmbSeedInd = 0
fgSeedInd = 1
phiSeedInd = 2
noiseSeedInd = 3
poissonSeedInd = 4

phiSeed = (0, 0, 2, i)
cmbSeed = (set, 0, 0, i)
fgseed = (set, 0, 1, i)

"""


class SimGen(object):
    def __init__(self, version, model="act_mr3", cmb_type='LensedCMB', dobeam=True, add_foregrounds=True, apply_window=True, max_cached=1,
                 extract_region = None,
                 extract_region_shape = None,
                 extract_region_wcs = None):
        
        """
        version: The version identifier for the filename of covsqrts on disk
        model: The name of an implemented soapack datamodel
        extract_region: An optional map whose footprint on to which the sims are extracted
        extract_region_shape: Instead of passing a map for extract_region, one can pass its shape and wcs
        extract_region_wcs: Instead of passing a map for extract_region, one can pass its shape and wcs
        max_cached: The maximum number of cached sim/or alms
        """
        self.noise_gen  = noise.NoiseGen(version=version,model=model,ncache=max_cached,verbose=False)
        self.signal_gen = signal.SignalGen(cmb_type=cmb_type, dobeam=dobeam, add_foregrounds=add_foregrounds, apply_window=apply_window, max_cached=max_cached, model=model)
        self.default_geometries = {}
        self.model = model
        
        if (extract_region is not None) or (extract_region_shape is not None):
            self._refoot = True
            if (extract_region is not None):
                assert extract_region_shape is None
                assert extract_region_wcs is None
                extract_region_shape, extract_region_wcs = extract_region.shape, extract_region.wcs
            assert extract_region_wcs is not None
            self._eshape, self._ewcs = extract_region_shape, extract_region_wcs
        else:
            self._refoot = False

    def get_default_geometry(self,season, patch, array, freq, mask_patch):
        if not patch in self.default_geometries.keys():
            oshape,owcs = self.noise_gen.load_covsqrt(season,patch,array,coadd=True,mask_patch=mask_patch,get_geometry=True)
            self.default_geometries[patch] = (oshape, owcs)
        else: pass
        return self.default_geometries[patch]
    
    def _footprint(self,imap):
        if not(self._refoot): return imap
        else: return enmap.extract(imap,self._eshape,self._wcs)

    ## Note: It will be probably better to use a decorator to delegate. For now, it does it explicetly

    def get_signal(self, season, patch, array, freq, sim_num, save_alm=True, save_map=False, set_idx=0,oshape=None,owcs=None,mask_patch=None,fgflux="15mjy", add_poisson_srcs=False):
        # return cmb+fg sim
        if oshape is None: oshape, owcs = self.get_default_geometry(season, patch, array, freq, mask_patch)
        return self._footprint(self.signal_gen.get_signal_sim(season, patch, array, freq, set_idx, sim_num, save_alm=save_alm, save_map=save_map,oshape=oshape,owcs=owcs, fgflux=fgflux, add_poisson_srcs=add_poisson_srcs))

    def get_cmb(self, season, patch, array, freq, sim_num, save_alm=False, set_idx=0,oshape=None,owcs=None,mask_patch=None):        
        if oshape is None: oshape, owcs = self.get_default_geometry(season, patch, array, freq, mask_patch) 
        return self._footprint(self.signal_gen.get_cmb_sim(season, patch, array, freq, set_idx, sim_num, save_alm=save_alm,oshape=oshape,owcs=owcs))

    def get_phi(self, season, patch, array, freq, sim_num, save_alm=False, set_idx=0,oshape=None,owcs=None,mask_patch=None): 
        if oshape is None: oshape, owcs = self.get_default_geometry(season, patch, array, freq, mask_patch)
        return self._footprint(self.signal_gen.get_phi_sim(patch, set_idx, sim_num, save_alm=save_alm, oshape=oshape,owcs=owcs))

    def get_kappa(self, season, patch, array, freq, sim_num, save_alm=False, set_idx=0,oshape=None,owcs=None,mask_patch=None): 
        if oshape is None: oshape, owcs = self.get_default_geometry(season, patch, array, freq, mask_patch)
        return self._footprint(self.signal_gen.get_kappa_sim(patch, set_idx, sim_num, save_alm=save_alm, oshape=oshape,owcs=owcs))
    
    def get_fg(self, season, patch, array, freq, sim_num, save_alm=False, set_idx=0,oshape=None,owcs=None,mask_patch=None,fgflux="15mjy"): 
        if oshape is None: oshape, owcs = self.get_default_geometry(season, patch, array, freq, mask_patch)
        return self._footprint(self.signal_gen.get_fg_sim(season, patch, array, freq, set_idx, sim_num, save_alm=save_alm, oshape=oshape,owcs=owcs, fgflux=fgflux))
    
    def get_noise(self, season=None,patch=None,array=None, sim_num=None,mask_patch=None,set_idx=0,apply_ivar=True):
        # indexing is slighly different for signal and noise sim code ..
        # array = None if (array is not None) and (freq is not None) else '{}_{}'.format(array, freq)   
        patch_id = None
        if mask_patch is not None: patch_id = int(mask_patch[-3:]) # handle s16 patches
        seed = (set_idx, 0, 3, sim_num) + self.noise_gen.dm.get_noise_sim_seed(season,patch,array,patch_id)
        return self._footprint(self.noise_gen.generate_sim(season,patch,array,seed=seed,mask_patch=mask_patch,apply_ivar=apply_ivar))

    def get_sim(self,season,patch,array,sim_num, save_alm=True, save_map=False, set_idx=0,mask_patch=None,fgflux="15mjy"):
        shape,wcs = self.noise_gen.load_covsqrt(season,patch,array,coadd=True,mask_patch=mask_patch,get_geometry=True)
        # (nfreqs,nsplits,npol,Ny,Nx)
        noises,ivars = self.get_noise(season=season,patch=patch,array=array, sim_num=sim_num,mask_patch=mask_patch,set_idx=set_idx,apply_ivar=False)
        #noises,ivars = self.get_noise(season=season,patch=patch,array=array, sim_num=0,mask_patch=mask_patch,set_idx=set_idx,apply_ivar=False)
        afreqs = self.noise_gen.dm.array_freqs[array]
        signals = []
        for afreq in afreqs:
            # This could be improved
            if self.model=='act_mr3': pfreq = afreq.split('_')[1] 
            elif self.model=='planck_hybrid': 
                pfreq = afreq
                assert season is None
                season = "planck"
                patch = "planck"
                array = "planck"

            imap = self.get_signal(season, patch, array, pfreq, sim_num, save_alm=save_alm, save_map=save_map, set_idx=set_idx,oshape=shape,owcs=wcs,fgflux=fgflux)
            signals.append(imap)
        owcs = imap.wcs
        # (nfreqs,npol,Ny,Nx)
        signals = enmap.enmap(np.stack(signals),owcs)[:,None,...]
        del imap
        assert wcsutils.equal(wcs,noises.wcs)
        assert wcsutils.equal(owcs,noises.wcs)
        return self._footprint(noise.apply_ivar_window(signals + noises,ivars))


def get_default_geometry(version, season, patch, array, freq, model='act_mr3'):
    noise_gen  = noise.NoiseGen(version=version,model=model,ncache=1,verbose=False)
    oshape,owcs = noise_gen.load_covsqrt(season,patch,array,coadd=True,mask_patch=None,get_geometry=True)
    return (oshape, owcs)

class Sehgal09Gen(SimGen):
    def __init__(self, version, model="act_mr3", cmb_type='LensedUnabberatedCMB', dobeam=True, add_foregrounds=True, apply_window=True, max_cached=1, extract_region = None,  extract_region_shape = None, extract_region_wcs = None, eulers=None):
        super(Sehgal09Gen, self).__init__(version=version, model=model, cmb_type=cmb_type, dobeam=dobeam, add_foregrounds=add_foregrounds,
                 apply_window=apply_window, max_cached=max_cached, 
                 extract_region = extract_region, extract_region_shape = extract_region_shape, extract_region_wcs = extract_region_wcs)

        self.signal_gen = signal.Sehgal09Gen(cmb_type=cmb_type, dobeam=dobeam, add_foregrounds=add_foregrounds, apply_window=apply_window, max_cached=max_cached, model=model,  eulers=eulers)

