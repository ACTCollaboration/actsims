import numpy as np, os
from pixell import enmap, fft, powspec, curvedsky
from . import flipperDict, simTools, act_data_model
import healpy as hp


input_param_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../inputParams')
input_path      = lambda x : os.path.join(input_param_dir, x)

signal_dict     = flipperDict.flipperDict()
noise_dict      = flipperDict.flipperDict()
signal_dict.read_from_file(input_path('signal.dict')) 
noise_dict.read_from_file(input_path('templateInputsMr3c.dict'))

cmb_dir         = os.path.join(os.path.dirname(os.path.abspath(__file__)), signal_dict['cmbDir'])
data_dir        = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../data')


class MR3PATCH_HELPER(object):
    # a helper class to quickly generate sims for given patch
    def __init__(self, set_idx, sim_idx, patch, cmb_type='LensedCMB', noisediag_only=False, dobeam=True, add_foregrounds=True, apply_window=True):
        self.patch    = patch.lower()
        assert(self.patch in ['boss', 'deep1', 'deep56', 'deep5', 'deep6', 'deep8'])
        assert(cmb_type in ['LensedCMB', 'UnlensedCMB', 'LensedUnabberatedCMB'])

        self.set_idx    = set_idx
        self.sim_idx    = sim_idx
        self.cmb_type   = cmb_type
        self.alms_base  = {}
        self.alms_patch = {}
        self.signals    = {}
        self.noises     = {}
        self.psas       = []
        self.psafs      = []
        freqs           = [] 
        self.template   = None
        self.dobeam     = dobeam
        for key in signal_dict['beamNames'].keys():
            if self.patch in key: 
                self.psafs.append(key)
                self.psas.append(key[:-5])
                freqs.append(key.split('_')[-1])
            else: pass
        self.psas  = list(set(self.psas))
        self.psas.sort()
        self.psafs.sort()
        self.freqs = list(set(freqs))
        self.freqs.sort()
        for psaf in self.psafs:
            self.noises[psaf] = {}
        self.noisediag_only  = noisediag_only
        self.add_foregrounds = add_foregrounds
        self.apply_window    = apply_window

    def get_signal_sim(self, psaf, save_alm=False, save_map=False):
        assert(psaf in self.psafs)
        patch, season, array, freq = psaf.split('_')
        print "loading sims for {}".format(psaf)
        if self.signals.has_key(psaf): 
            print "loading precomputed sim {}".format(psaf)
            return self.signals[psaf].copy()
        
        alm_patch = None
        if self.alms_patch.has_key(psaf):  
            print "loading precomputed alm {}".format(psaf)
            return self.signals[psaf].copy()
            alm_patch = self.alms_patch[psaf].copy()
        else:
            if len(self.alms_base) == 0: self.load_alms_base()
            alm_patch = self.alms_base[freq].copy()
            if self.dobeam:
                print "apply beam for alm {}".format(psaf)
                _, beam_data = act_data_model.load_normalized_beam(season, array, patch, freq)
                for idx in range(alm_patch.shape[0]):
                    alm_patch[idx] = hp.sphtfunc.almxfl(alm_patch[idx].copy(), beam_data)
            else: pass
            if save_alm: self.alms_patch[psaf] = alm_patch.copy()

        signal = self.get_template()
        curvedsky.alm2map(alm_patch, signal, spin = [0,2], verbose=True)
        
        if self.apply_window:
            print 'apply window'
            axes = [-2, -1]
            for idx in range(signal.shape[0]):
                kmap   = fft.fft(signal[idx], axes=axes)
                wy, wx = enmap.calc_window(kmap.shape)
                wind   = wy[:,None]**1 * wx[None,:]**1
                kmap   *= wind

                signal[idx] = (fft.ifft(kmap, axes=axes, normalize=True)).real
                del kmap

        if save_map: self.signals[psaf] = signal.copy()
        return signal

    def get_noise_sim(self, psaf, seed=None, save_map=False):
        assert(psaf in self.psafs)
        ## should we preload noise templates

        if seed is None: seed = self.sim_idx
        patch, season, array, freq = psaf.split('_')
        print "loading noise sims for {} seed {}".format(psaf, seed) 

        if seed is not None and self.noises[psaf].has_key(seed): 
            print "loading precomputed sim for {} seed {}".format(psaf, seed)
            return self.noises[psaf][seed].copy()
        else:
            noise = simTools.getActpolSim(iterationNum=seed,
                        patch=patch, season=season, array=array, simType='noise', 
                        noiseDiagsOnly=self.noisediag_only, cmbSet=self.set_idx)
            if save_map:
                if 'pa3' not in psaf:
                    self.noises[psaf][seed] = noise[0].copy()
                else:
                    psaf_temp = '%s_{}' %('_'.join([patch,season,array]))
                    self.noises[psaf_temp.format('f090')][seed] = noise[0].copy()
                    self.noises[psaf_temp.format('f150')][seed] = noise[1].copy()

            freq_idx = 1 if noise.shape[0] == 2 and freq == 'f150' else 0
            noise = noise[freq_idx]
            return noise
 

    def get_template(self):
        if self.template is None:
            template_file = os.path.join(os.path.join(os.path.dirname(os.path.abspath(__file__))), noise_dict['dataMapDir'], 
                    'totalWeightMapI_{}_fromenlib.fits'.format(self.psafs[0]))
            template      = enmap.read_map(template_file)
            self.template = enmap.empty((3,) + template.shape, template.wcs)
        else: pass
        return self.template.copy()

    def load_alms_base(self):
        # note: beam is set to false
        print("loading alm base")
        cmb_file   = os.path.join(cmb_dir, 'fullsky%s_alm_set%02d_%05d.fits' %(self.cmb_type, self.set_idx, self.sim_idx))
        print("loading %s" %cmb_file)
        alm_signal = np.complex128(hp.fitsfunc.read_alm(cmb_file, hdu = (1,2,3))) 
        alm_signal = np.tile(alm_signal, (len(self.freqs), 1, 1))
        if self.add_foregrounds:
            print("adding fgs to the base")
            seed         = (self.set_idx, 0, 1, self.sim_idx)# copying the structure from simtools
            fg_file      = os.path.join(data_dir, signal_dict['foregroundPowerFile'])
            fg_power     = powspec.read_spectrum(fg_file, ncol = 3, expand = 'row')
            
            lmax_sg      = hp.Alm.getlmax(alm_signal.shape[-1])
            alm_fg90_150 = curvedsky.rand_alm(fg_power, seed = seed)#, lmax=lmax_sg)
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

        for i, freq in enumerate(self.freqs):
            self.alms_base[freq] = alm_signal[i].copy()
        del alm_signal

    def clear(self):
        for key in self.alms_base.keys():
            del self.alms_base[key]
        for key in self.alms_patch.keys():
            del self.alms_patch[key]
        for key in self.signals.keys():
            del self.signals[key]
        for key in self.noises.keys():
            del self.noises[key]


