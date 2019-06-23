import numpy as np, os
from pixell import enmap, powspec, curvedsky, fft as pfft
from soapack import interfaces as sints
import healpy as hp
import warnings
from collections import OrderedDict as ODict
from actsims.util import SEED_TRACKER as seedgen
from itertools import combinations
actsim_root = os.path.dirname(os.path.realpath(__file__))

class SignalGen(object):
    # a helper class to quickly generate sims for given patch
    def __init__(self, cmb_type='LensedCMB', dobeam=True, add_foregrounds=True, apply_window=True, max_cached=1, model="act_mr3", apply_rotation=False, alpha_map=None):
        """
        model: The name of an implemented soapack datamodel
        ncache: The number of 

        """ 
        self.data_model = sints.models[model]()
        self.cmb_types   = ['LensedCMB', 'UnlensedCMB', 'LensedUnabberatedCMB']
        paths            = sints.dconfig['actsims']
        self.signal_path = paths['signal_path']
        assert(self.signal_path is not None)
        assert(cmb_type in self.cmb_types)

        self.supported_sims = []
        # ACT
        patches = {'s13':['deep1','deep5','deep6'],'s14':['deep56'],'s15':['deep56','boss','deep8'],'s16':['cmb']}
        arrays = {'s13':['pa1_f150'],'s14':['pa1_f150','pa2_f150'],'s15':['pa1_f150','pa2_f150','pa3_f150','pa3_f090'],'s16':['pa2_f150','pa3_f150','pa3_f090']}
        for season in patches.keys():
            for patch in patches[season]:
                for array in arrays[season]:
                    sstr = "%s_%s_%s" % (season,patch,array)
                    self.supported_sims.append(sstr)

        # Planck
        for freq in [30,44,70,100,143,217,353,545,857]:
            self.supported_sims.append("planck_planck_planck_%s" % (str(freq).zfill(3)))

        self.supported_sims.sort()
        self.freqs           = ['f150','f090']
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
        self.apply_rotation   = apply_rotation
        self.alpha_map        = alpha_map

    def is_supported(self, season, patch,array, freq):
        signal_idx      = self.__combine_idxes__(season, patch, array , freq)
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

    def get_signal_sim(self, season, patch, array, freq, set_idx, sim_num, oshape, owcs, save_alm=False, save_map=False, fgflux="15mjy", add_poisson_srcs=False):
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
                self.load_alms_base(set_idx, sim_num, fgflux=fgflux)
            alm_patch = self.alms_base[base_alm_idx][freq_idx].copy()
            if add_poisson_srcs: 
                alm_patch[0] += self.get_poisson_srcs_alms(set_idx, sim_num, patch, alm_patch[0].shape, oshape=oshape, owcs=owcs)
            if self.dobeam:
                print ("apply beam for alm {}".format(signal_idx))
                alm_patch = self.__apply_beam__(alm_patch, season, patch, array, freq)
            else: pass
            if save_alm: 
                self.manage_cache(self.alms_patch, self.max_cached-1) 
                self.alms_patch[signal_idx] = alm_patch.copy()
        
        return self.__signal_postprocessing__(patch, signal_idx, alm_patch, save_map=save_map,oshape=oshape,owcs=owcs, apply_window=self.apply_window)
 
    def get_cmb_sim(self, season, patch, array, freq, set_idx, sim_num, oshape, owcs, save_alm=False ):
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


    def get_fg_sim(self, season, patch, array, freq, set_idx, sim_num, oshape, owcs, save_alm=False, fgflux="15mjy"):
        assert(self.is_supported(season, patch, array, freq))
        print("[WARNING] get_fg_sim() is implemented for debugging purpose. Use get_signal_sim() for the production run")

        signal_idx   = self.get_signal_idx(season, patch, array, freq, set_idx, sim_num)
        alm_fg = None
        if signal_idx in self.alms_fg:  
            print ("loading precomputed alm cmb {}".format(signal_idx))
            alm_fg = self.alms_fg[signal_idx].copy()
        else:
            alm_fg90_150  = self.load_alm_fg(set_idx, sim_num, fgflux=fgflux) 
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
        
    def get_phi_sim(self, patch, set_idx, sim_num, oshape, owcs, save_alm=False):
        return self.__get_lens_potential_sim__(patch, set_idx, sim_num, mode='phi', save_alm=save_alm, oshape=oshape, owcs=owcs)

    def get_kappa_sim(self, patch, set_idx, sim_num, oshape, owcs, save_alm=False): 
        return self.__get_lens_potential_sim__(patch, set_idx, sim_num, mode='kappa', save_alm=save_alm, oshape=oshape, owcs=owcs)
    
    def __get_lens_potential_sim__(self, patch, set_idx, sim_num, oshape, owcs, mode='phi', save_alm=False):
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
        beam_data = self.data_model.get_beam(l_beam, season=season, patch=patch, array='{}_{}'.format(array, freq) if array!='planck' else freq) 
        
        for idx in range(alm_patch.shape[0]):
            alm_patch[idx] = hp.sphtfunc.almxfl(alm_patch[idx].copy(), beam_data)
        return alm_patch 

    def __signal_postprocessing__(self, patch, signal_idx, alm_patch, save_map, oshape, owcs, apply_window=True):
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


    def load_alms_base(self, set_idx, sim_idx, cache=True, fg_override=None, ret_alm=False, fgflux="15mjy"):
        # note: beam is set to false
        print("loading alm base")
        cmb_file   = os.path.join(self.signal_path, 'fullsky%s_alm_set%02d_%05d.fits' %(self.cmb_type, set_idx, sim_idx))
        print("loading %s" %cmb_file)
        alm_signal = np.complex128(hp.fitsfunc.read_alm(cmb_file, hdu = (1,2,3))) 

        # check if rotation is needed
        if self.apply_rotation:
            if not np.any(self.alpha_map):
                print("[Warning] want to apply rotation but alpha map is not provided...")
            else:
                alpha_map = self.alpha_map
                print("applying rotation field")
                # get wcs from rot_map: doesn't matter which wcs is used because eventually
                # it is converted back to alm
                wcs = alpha_map.wcs
                shape = (3,) + alpha_map.shape[-2:]

                # generate cmb map from alm_signal
                cmb_map = enmap.empty(shape, wcs)
                curvedsky.alm2map(alm_signal, cmb_map)

                # apply the rotation field
                cmb_map = enmap.rotate_pol(cmb_map, 2*alpha_map)

                # convert back to alm
                lmax = hp.Alm.getlmax(alm_signal.shape[-1])
                alm_signal = curvedsky.map2alm(cmb_map, lmax=lmax)
                del cmb_map

        alm_signal = np.tile(alm_signal, (len(self.freqs), 1, 1))

        if fg_override is not None: print("[Warning] override the default foreground flag....")
        add_foregrounds = fg_override if fg_override is not None else self.add_foregrounds
        if add_foregrounds:
            print("adding fgs to the base")     
            alm_fg90_150 = self.load_alm_fg(set_idx, sim_idx, fgflux=fgflux)
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
    
    def load_alm_fg(self, set_idx, sim_idx, fgflux):
        print("loading fg alm")
        
        seed = seedgen.get_fg_seed(set_idx, sim_idx, fgflux)
        if fgflux == "15mjy":
            print("loading FG with 15mJy fluxcut")
            fg_file      = os.path.join(actsim_root, '../data/fg.dat')
        elif fgflux=="100mjy":
            print("loading FG with 100mJy fluxcut")
            fg_file      = os.path.join(actsim_root, '../data/highflux_fg.dat')
        else:
            assert(False) ### :o
        fg_power     = powspec.read_spectrum(fg_file, ncol = 3, expand = 'row')
        alm_fg90_150 = curvedsky.rand_alm_healpy(fg_power, seed = seed)#, lmax=lmax_sg)
        return alm_fg90_150
    
    def get_template(self, patch, shape, wcs):
        if patch not in self.templates:
            self.manage_cache(self.templates, self.max_cached-1) 
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


    def get_poisson_srcs_alms(self, set_idx, sim_num, patch, alm_shape, oshape, owcs):

        def deltaTOverTcmbToJyPerSr(freqGHz,T0 = 2.726):
            """
            @brief the function name is self-eplanatory
            @return the converstion factor
            stolen from Flipper -- van engelen
            """
            kB = 1.380658e-16
            h = 6.6260755e-27
            c = 29979245800.
            nu = freqGHz*1.e9
            x = h*nu/(kB*T0)
            cNu = 2*(kB*T0)**3/(h**2*c**2)*x**4/(4*(np.sinh(x/2.))**2)
            cNu *= 1e23
            return cNu

        TCMB_uk = 2.72e6
        
        if oshape[0] > 3:
            #then this is a multichroic array, and sadly we only have this at 150 GHz for now
            raise Exception('get_poisson_srcs_alms only implemented for 150 GHz so far ' \
                            + '(that is the model we currently have for radio sources) ')
        else:
            freq_ghz = 148
        
        #ideally this RNG stuff would be defined in a central place to
        #avoid RNG collisions.  Old version is currently commented out at top of
        #simgen.py
        templ = self.get_template(patch, shape = oshape, wcs = owcs)
        
        templ[:] = 0
        seed = seedgen.get_poisson_seed(set_idx, sim_num)
        np.random.seed(seed = seed)

        #Wasn't sure how to codify this stuff outside this routine - hardcoded for now
        S_min_Jy = .001
        S_max_Jy = .015


        tucci = np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../data/ns_148GHz_modC2Ex.dat'))

        S = tucci[:, 0]
        dS = S[1:] - S[0:-1]
        dS = np.append(dS, [0.])
        dNdS = tucci[:, 1]

        mean_numbers_per_patch = dNdS * enmap.area(templ.shape, templ.wcs) * dS

        numbers_per_fluxbin = np.random.poisson(mean_numbers_per_patch)

        #note pixel areas not constant for pixell maps
        pixel_areas = enmap.pixsizemap(templ.shape, templ.wcs)
        
        for si, fluxval in enumerate(S[S <= S_max_Jy]):
            xlocs = np.random.randint( templ.shape[-1], size = numbers_per_fluxbin[si])
            ylocs = np.random.randint( templ.shape[-2], size = numbers_per_fluxbin[si])

            #add the value in jy / sr, i.e. divide by the solid angle of a pixel.
            templ[0, ylocs, xlocs] += fluxval / pixel_areas[ylocs, xlocs]

        map_factor = TCMB_uk / deltaTOverTcmbToJyPerSr(freq_ghz)
        templ *= map_factor

        #GET ALMs
        import healpy
        output = curvedsky.map2alm(templ[0], lmax = healpy.Alm.getlmax(alm_shape[0]))
        # import pdb
        # pdb.set_trace()

        return output
    
class SRCFREE_SPECS(object):
    def __init__(self):
        self.supported_psa  = {'deep56': [], 'boss': []}
        self.supported_spec = {'deep56': [], 'boss': []}
        self.act_arrays     = {'deep56':['d56_%s' % str(x).zfill(2) for x in range(1,7) ], \
        'boss':['boss_%s'%str(x).zfill(2) for x in range(1,4) ]}
        self.planck_arrays  = ['p0%d' % x for x in range(1,9)]
        self.full_arrays    = {'deep56': [], 'boss': []}
        for patch in self.supported_psa.keys():
            self.full_arrays[patch] =  self.act_arrays[patch] + self.planck_arrays
        self.data_models    = {'act_mr3':'act','planck_hybrid':'planck'}
        self.spec_storage   = None
        for patch in self.supported_psa.keys():
            for array in self.full_arrays[patch]:
                dm        = self.data_models[sints.arrays(array,'data_model')]
                season    = sints.arrays(array,'season')
                array_str = '{}_{}'.format(sints.arrays(array,'array'),sints.arrays(array,'freq'))
                if dm is 'planck': 
                    array_str = sints.arrays(array,'freq')
                    # filter out planck low-fi (below 100GHz)
                    is_lowfi = lambda x: float(x) < 100
                    if is_lowfi(array_str): continue
                self.supported_psa[patch].append(self.__single_idx__(dm, season, array_str))
            self.supported_psa[patch].sort()

        # ordering (1) actxact, (2) actxplanck, (3) planckxplanck 
        def __append_supported_spec_idx__(exp1, exp2):
            for patch in self.supported_psa.keys():
                for i, idxi in enumerate(self.supported_psa[patch]):
                    for j, idxj in enumerate(self.supported_psa[patch][i:]):
                        if not(exp1 in idxi and exp2 in idxj): continue
                        spec_idx = '{}_{}_{}'.format(patch, idxi, idxj)
                        self.supported_spec[patch].append(spec_idx)
        __append_supported_spec_idx__('act','act')
        __append_supported_spec_idx__('act','planck')
        __append_supported_spec_idx__('plack','planck')

        data_file          = os.path.join(actsim_root, '../data/spectra/srcfre_fit_mr3f_cal.npz')
        self.spec_storage  = np.load(data_file)
        self.covs          = {}

    def __single_idx__(self, dm, season, array):
        return '{}_{}_{}'.format(dm, season, array)

    def __double_idx__(self, dm1, season1, array1, dm2, season2, array2):
        idx1 = self.__single_idx__(dm1, season1, array1)
        idx2 = self.__single_idx__(dm2, season2, array2)

        temp = [idx1,idx2]
        temp.sort()
        return '{}_{}'.format(temp[0], temp[1])

    def get_spec_idx(self, patch, dm1, season1, array1, dm2, season2, array2):
        return '{}_{}'.format(patch, self.__double_idx__(dm1, season1, array1, dm2, season2, array2))
    def is_supported(self, patch, dm1, season1, array1, dm2, season2, array2):
        spec_idx = self.get_spec_idx(patch, dm1, season1, array1, dm2, season2, array2)
        return spec_idx in  self.supported_spec[patch]

    def get_cov(self, patch):
        if patch in self.covs.keys(): return self.covs[patch].copy() # return cached cov
        npad  = int(self.spec_storage['l'][0]) # should have started at ell=0, but it starts at ell=2
        nelmt = len(self.spec_storage['l'])+npad
        nspec = len(self.supported_psa[patch])
        cov   = np.zeros((nspec, nspec, nelmt))
        
        temp = lambda x, y, z: '{}_{}_{}'.format(x,y,z)
        for i, psai in enumerate(self.supported_psa[patch]):
            for j, psaj in enumerate(self.supported_psa[patch][i:]):
                cov[i,j][npad:] = self.spec_storage[temp(patch, psai, psaj)].copy()
                if i != j: cov[j,i][npad:] = self.spec_storage[temp(patch, psai, psaj)].copy()
        self.covs[patch] = cov
        return self.covs[patch]


        





