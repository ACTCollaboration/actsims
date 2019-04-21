import numpy as np, os
from pixell import enmap, powspec, curvedsky, fft as pfft
from soapack import interfaces as sints
import healpy as hp
import warnings
from collections import OrderedDict as ODict

actsim_root = os.path.dirname(os.path.realpath(__file__))

class SignalGen(object):
    # a helper class to quickly generate sims for given patch
    def __init__(self, cmb_type='LensedCMB', dobeam=True, add_foregrounds=True, apply_window=True, max_cached=1, model="act_mr3", extract_region=None,extract_region_shape=None, extract_region_wcs=None):
        """
        model: The name of an implemented soapack datamodel
        extract_region: An optional map whose footprint on to which the sims are made
        extract_region_shape: Instead of passing a map for extract_region, one can pass its shape and wcs
        extract_region_wcs: Instead of passing a map for extract_region, one can pass its shape and wcs
        ncache: The number of 

        """ 
        self.data_model = sints.models[model](region=extract_region, region_shape=extract_region_shape,region_wcs=extract_region_wcs)
        self.cmb_types   = ['LensedCMB', 'UnlensedCMB', 'LensedUnabberatedCMB']
        paths            = sints.dconfig['actsims']
        self.signal_path = paths['signal_path']
        assert(self.signal_path is not None)
        assert(cmb_type in self.cmb_types)


        self.supported_sims = ['s13_pa1_deep1_f150', 's13_pa1_deep5_f150', 's13_pa1_deep6_f150', 's14_pa1_deep56_f150', 's14_pa2_deep56_f150', \
                               's15_pa1_boss_f150', 's15_pa1_deep56_f150', 's15_pa1_deep8_f150', 's15_pa2_boss_f150', 's15_pa2_deep56_f150', \
                               's15_pa2_deep8_f150','s15_pa3_boss_f090', 's15_pa3_boss_f150', 's15_pa3_deep56_f090', 's15_pa3_deep56_f150', \
                               's15_pa3_deep8_f090', 's15_pa3_deep8_f150','s16_pa2_cmb_f150','s16_pa3_cmb_f150','s16_pa3_cmb_f090']
        self.supported_sims.sort()
        self.freqs           = ['f090', 'f150']
        self.cmb_type         = cmb_type
        self.max_cached       = max_cached
        self.alms_base        = ODict()
        self.alms_cmb         = ODict()
        self.alms_phi         = ODict()
        self.alms_kappa       = ODict()
        self.alms_fg          = ODict()
        self.alms_patch       = ODict()
        self.signals          = ODict()
        self.templates        = ODict()
        self.apply_window     = apply_window
        self.add_foregrounds  = add_foregrounds
        self.dobeam           = dobeam

    def is_supported(self, sesaon, array, patch, freq):
        signal_idx      = self.__combine_idxes__(sesaon, patch, array , freq)
        supported = signal_idx in self.supported_sims
        if not supported:
            print("unknown type of sims: {} ".format(signal_idx))
            print("supported sims are {}".format(self.supported_sims))

        return supported

    def __combine_idxes__(self, season, patch, array, freq):
        return '_'.join([season, patch, array, freq])
    
    def get_base_alm_idx(self, set_idx, sim_num):
        return '_'.join(['set0%d'%set_idx, '%05d'%sim_num])

    def get_signal_idx(self, season, patch, array, freq, set_idx, sim_num):
        return '_'.join([season, patch, array, freq, 'set0%d'%set_idx, '%05d'%sim_num])

    def get_signal_sim(self, season, patch, array, freq, set_idx, sim_num, save_alm=False, save_map=False,oshape=None,owcs=None):
        assert(self.is_supported(season, patch, array, freq))

        base_alm_idx = self.get_base_alm_idx(set_idx, sim_num) 
        signal_idx   = self.get_signal_idx(season, patch, array, freq, set_idx, sim_num)

        print( "loading sims for {}".format(signal_idx))
        if signal_idx in self.signals: 
            print ("loading precomputed sim {}".format(signal_idx))
            return self.signals[signal_idx].copy()
        
        alm_patch = None
        if signal_idx in self.alms_patch:  
            print ("loading precomputed alm {}".format(signal_idx))
            alm_patch = self.alms_patch[signal_idx].copy()
        else:
            freq_idx  = 0 if freq == 'f090' else 1
            if base_alm_idx not in self.alms_base:
                self.manage_cache(self.alms_base, self.max_cached-1)
                self.load_alms_base(set_idx, sim_num)
            alm_patch = self.alms_base[base_alm_idx][freq_idx].copy()
            if self.dobeam:
                print ("apply beam for alm {}".format(signal_idx))
                alm_patch = self.__apply_beam__(alm_patch, season, patch, array, freq)
            else: pass
            if save_alm: 
                self.manage_cache(self.alms_patch, self.max_cached-1) 
                self.alms_patch[signal_idx] = alm_patch.copy()
        
        return self.__signal_postprocessing__(patch, signal_idx, alm_patch, save_map=save_map,oshape=oshape,owcs=owcs, apply_window=self.apply_window)
 
    def get_cmb_sim(self, season, patch, array, freq, set_idx, sim_num, save_alm=False,oshape=None,owcs=None):
        assert(self.is_supported(season, patch, array, freq))
        print("[WARNING] get_cmb_sim() is implemented for debugging purpose. Use get_signal_sim() for the production run")

        signal_idx   = self.get_signal_idx(season, patch, array, freq, set_idx, sim_num)
        
        alm_cmb = None
        if signal_idx in self.alms_cmb:  
            print ("loading precomputed alm cmb {}".format(signal_idx))
            alm_cmb = self.alms_cmb[signal_idx].copy()
        else:
            alm_cmb = self.load_alms_base(set_idx, sim_num, cache=False, fg_override=False, ret_alm=True)
            freq_idx  = 0 if freq == 'f090' else 1
            alm_cmb = alm_cmb[freq_idx].copy()
            if self.dobeam:
                print ("apply beam for alm {}".format(signal_idx))
                alm_cmb = self.__apply_beam__(alm_cmb, season, patch, array, freq)
            else: pass
            if save_alm: 
                self.manage_cache(self.alms_cmb, self.max_cached-1) 
                self.alms_cmb[signal_idx] = alm_cmb.copy()

        return self.__signal_postprocessing__(patch, signal_idx, alm_cmb, save_map=False,oshape=oshape,owcs=owcs, apply_window=self.apply_window)


    def get_fg_sim(self, season, patch, array, freq, set_idx, sim_num, save_alm=False,oshape=None,owcs=None):
        assert(self.is_supported(season, patch, array, freq))
        print("[WARNING] get_fg_sim() is implemented for debugging purpose. Use get_signal_sim() for the production run")

        signal_idx   = self.get_signal_idx(season, patch, array, freq, set_idx, sim_num)
        alm_fg = None
        if signal_idx in self.alms_fg:  
            print ("loading precomputed alm cmb {}".format(signal_idx))
            alm_fg = self.alms_fg[signal_idx].copy()
        else:
            alm_fg90_150  = self.load_alm_fg(set_idx, sim_num) 
            alm_out = np.zeros((3, len(alm_fg90_150[-1])), dtype = np.complex128) 

            freq_idx      = 1 if freq == 'f150' else 0
            alm_out[0, :] = alm_fg90_150[freq_idx, :].copy()
       
            alm_fg = alm_out
            if self.dobeam:
                print ("apply beam for alm {}".format(signal_idx))
                alm_fg = self.__apply_beam__(alm_fg, season, patch, array, freq)
            else: pass
            if save_alm: 
                self.manage_cache(self.alms_fg, self.max_cached-1) 
                self.alms_fg[signal_idx] = alm_fg.copy()

        return self.__signal_postprocessing__(patch, signal_idx, alm_fg, save_map=False,oshape=oshape,owcs=owcs, apply_window=self.apply_window)
        
    def get_phi_sim(self, patch, set_idx, sim_num, save_alm=False, oshape=None, owcs=None):
        return self.__get_lens_potential_sim__(patch, set_idx, sim_num, mode='phi', save_alm=save_alm, oshape=oshape, owcs=owcs)

    def get_kappa_sim(self, patch, set_idx, sim_num, save_alm=False, oshape=None, owcs=None): 
        return self.__get_lens_potential_sim__(patch, set_idx, sim_num, mode='kappa', save_alm=save_alm, oshape=oshape, owcs=owcs)
    
    def __get_lens_potential_sim__(self, patch, set_idx, sim_num, mode='phi', save_alm=False, oshape=None, owcs=None):
        assert(mode in ['phi', 'kappa'])
        lenp_idx     = self.get_base_alm_idx(set_idx, sim_num)
        alms_lenp    = self.alms_phi     if mode == 'phi' else self.alms_kappa       
        loader_func  = self.load_alm_phi if mode == 'phi' else self.load_alm_kappa 

        alm_lenp = None
        if lenp_idx in alms_lenp:  
            print ("loading precomputed alm lenp {}".format(lenp_idx))
            alm_lenp = alms_lenp[alm_lenp].copy()
        else:
            alm_lenp = loader_func(set_idx, sim_num, cache=False, ret_alm=True)
            if save_alm: 
                self.manage_cache(alms_lenp, self.max_cached-1) 
                alms_lenp[lenp_idx] = alm_lenp.copy()

        return self.__signal_postprocessing__(patch, lenp_idx, alm_lenp, save_map=False, oshape=oshape,owcs=owcs, apply_window=False)

    def __apply_beam__(self, alm_patch, season, patch, array, freq):
        lmax      = hp.Alm.getlmax(alm_patch.shape[-1])
        l_beam    = np.arange(0, lmax+100, dtype=np.float)
        beam_data = self.data_model.get_beam(l_beam, season, patch, '{}_{}'.format(array, freq)) 
        
        for idx in range(alm_patch.shape[0]):
            alm_patch[idx] = hp.sphtfunc.almxfl(alm_patch[idx].copy(), beam_data)
        return alm_patch 

    def __signal_postprocessing__(self, patch, signal_idx, alm_patch, save_map, oshape=None, owcs=None, apply_window=True):
        signal = self.get_template(patch,shape=oshape,wcs=owcs)
        signal = signal if len(alm_patch.shape) > 1 else signal[0,...]
        curvedsky.alm2map(alm_patch, signal, spin = [0,2], verbose=True)

        if apply_window:
            print('apply window')
            axes = [-2, -1]
            for idx in range(signal.shape[0]):
                kmap   = pfft.fft(signal[idx], axes=axes)
                wy, wx = enmap.calc_window(kmap.shape)
                wind   = wy[:,None]**1 * wx[None,:]**1
                kmap   *= wind

                signal[idx] = (pfft.ifft(kmap, axes=axes, normalize=True)).real
                del kmap

        if save_map: 
            self.manage_cache(self.signals, self.max_cached-1)
            self.signals[signal_idx] = signal.copy()
        return signal


    def load_alms_base(self, set_idx, sim_idx, cache=True, fg_override=None, ret_alm=False):
        # note: beam is set to false
        print("loading alm base")
        cmb_file   = os.path.join(self.signal_path, 'fullsky%s_alm_set%02d_%05d.fits' %(self.cmb_type, set_idx, sim_idx))
        print("loading %s" %cmb_file)
        alm_signal = np.complex128(hp.fitsfunc.read_alm(cmb_file, hdu = (1,2,3))) 
        alm_signal = np.tile(alm_signal, (len(self.freqs), 1, 1))

        if fg_override is not None: print("[Warning] override the default foreground flag....")
        add_foregrounds = fg_override if fg_override is not None else self.add_foregrounds
        if add_foregrounds:
            print("adding fgs to the base")     
            alm_fg90_150 = self.load_alm_fg(set_idx, sim_idx)
            lmax_sg      = hp.Alm.getlmax(alm_signal.shape[-1])
            alm_out      = np.zeros((len(self.freqs), 3, len(alm_fg90_150[-1])), dtype = np.complex128) 
            lmax_fg      = hp.Alm.getlmax(alm_fg90_150.shape[-1])

            for idx, freq in enumerate(self.freqs):
                freq_idx = 1 if freq == 'f150' else 0
                alm_out[idx, 0, :] = alm_fg90_150[freq_idx, :].copy()
            
            for m in range(lmax_sg+1):
                lmin = m
                lmax = lmax_sg

                idx_ssidx = hp.Alm.getidx(lmax_sg, lmin, m)
                idx_seidx = hp.Alm.getidx(lmax_sg, lmax, m)
                idx_fsidx = hp.Alm.getidx(lmax_fg, lmin, m)
                idx_feidx = hp.Alm.getidx(lmax_fg, lmax, m)

                alm_out[..., idx_fsidx:idx_feidx+1] = alm_out[..., idx_fsidx:idx_feidx+1] + alm_signal[..., idx_ssidx:idx_seidx+1]

            alm_signal = alm_out.copy()
            del alm_out, alm_fg90_150

        if cache: self.alms_base[self.get_base_alm_idx(set_idx, sim_idx)] = alm_signal.copy()
        
        if ret_alm: return alm_signal
        del alm_signal

    def load_alm_phi(self, set_idx, sim_idx, cache=True, ret_alm=False):
        # note: beam is set to false
        print("loading alm phi")
        phi_file = os.path.join(self.signal_path, 'fullskyPhi_alm_set%02d_%05d.fits' %(set_idx, sim_idx))
        print("loading %s" %phi_file)
        alm_phi  = np.complex128(hp.fitsfunc.read_alm(phi_file, hdu = (1))) 

        if cache: self.alms_phi[self.get_base_alm_idx(set_idx, sim_idx)] = alm_phi.copy()
        
        if ret_alm: return alm_phi
        del alm_phi
    
    def load_alm_kappa(self, set_idx, sim_idx, cache=True, ret_alm=False):
        # note: beam is set to false
        alm_phi  = self.load_alm_phi(set_idx, sim_idx, False, True)
        print("alm_phi to alm_kappa")
        lmax = hp.Alm.getlmax(alm_phi.shape[-1])
        l    = np.arange(0, lmax+1, dtype=np.float)

        alm_kappa = hp.sphtfunc.almxfl(alm_phi, l*(l+1.)/2.)
        if cache: self.alms_kappa[self.get_base_alm_idx(set_idx, sim_idx)] = alm_kappa.copy()
        
        if ret_alm: return alm_kappa
        del alm_kappa
    
    def load_alm_fg(self, set_idx, sim_idx):
        print("loading fg alm")
        seed         = (set_idx, 0, 1, sim_idx)# copying the structure from simtools
        fg_file      = os.path.join(actsim_root, '../data/fg.dat')
        fg_power     = powspec.read_spectrum(fg_file, ncol = 3, expand = 'row')
        print(fg_power.shape,seed)
        alm_fg90_150 = curvedsky.rand_alm_healpy(fg_power, seed = seed)#, lmax=lmax_sg)
        return alm_fg90_150
    
    def get_template(self, patch,shape=None,wcs=None):
        if patch not in self.templates:
            self.manage_cache(self.templates, self.max_cached-1) 
            if shape is None:
                template      = self.data_model.get_mask(patch)
                shape,wcs = template.shape,template.wcs
            self.templates[patch] = enmap.empty((3,) + shape[-2:], wcs)
        else: pass
        return self.templates[patch].copy()

    def manage_cache(self, odict, max_cached=None):
        if max_cached is None: max_cached = self.max_cached
        if max_cached < 0: max_cached = 0
        nelmt = len(odict)
        for key in odict.keys():
            if nelmt <= max_cached: continue
            print("purging {} from the cache".format(key))
            del odict[key]; nelmt -= 1

    def clear(self):
        for key in self.alms_base.keys():
            del self.alms_base[key]
        for key in self.alms_patch.keys():
            del self.alms_patch[key]
        for key in self.signals.keys():
            del self.signals[key]


