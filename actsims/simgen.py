from actsims import signal, noise

class SimGen(object):
    def __init__(self, version, model="act_mr3", extract_region=None,extract_region_shape=None, extract_region_wcs=None, cmb_type='LensedCMB', dobeam=True, add_foregrounds=True, apply_window=True, max_cached=1):
        
        """
        version: The version identifier for the filename of covsqrts on disk
        model: The name of an implemented soapack datamodel
        extract_region: An optional map whose footprint on to which the sims are made
        extract_region_shape: Instead of passing a map for extract_region, one can pass its shape and wcs
        extract_region_wcs: Instead of passing a map for extract_region, one can pass its shape and wcs
        max_cached: The maximum number of cached sim/or alms
        """
 
        
        self.noise_gen  = noise.NoiseGen(version=version,model=model,extract_region=extract_region,extract_region_shape=extract_region_shape,extract_region_wcs=extract_region_wcs,ncache=max_cached,verbose=False)
        self.signal_gen = signal.SignalGen(cmb_type=cmb_type, dobeam=dobeam, add_foregrounds=add_foregrounds, apply_window=apply_window, max_cached=max_cached, model=model, extract_region=extract_region,extract_region_shape=extract_region_shape, extract_region_wcs=extract_region_wcs)


    ## Note: It will be probably better to use a decorator to delegate. For now, it does it explicetly

    def get_signal(self, season, array, patch, freq, set_idx, sim_num, save_alm=True, save_map=False):
        # return cmb+fg sim
        return self.signal_gen.get_signal_sim(season, array, patch, freq, set_idx, sim_num, save_alm, save_map)

    def get_cmb(self, season, array, patch, freq, set_idx, sim_num, save_alm=False):
        return self.signal_gen.get_cmb_sim(season, array, patch, freq, set_idx, sim_num, save_alm)

    def get_fg(self, season, array, patch, freq, set_idx, sim_num, save_alm=False):
        return self.signal_gen.get_fg_sim(season, array, patch, freq, set_idx, sim_num, save_alm)

    def get_noise(self, season=None,patch=None,array=None, freq=None, seed=None,mask_patch=None,binary_percentile=10.):
        # indexing is slighly different for signal and noise sim code ..
        array = None if (array is not None) and (freq is not None) else '{}_{}'.format(array, freq)   
        return self.noise_gen.generate_sim(season,patch,array,seed,mask_patch,binary_percentile)

