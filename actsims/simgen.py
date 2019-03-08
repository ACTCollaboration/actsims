from actsims import signal, noise

class SimGen(object):
    def __init__(self, version, model="act_mr3", extract_region=None,extract_region_shape=None, extract_region_wcs=None, cmb_type='LensedCMB', dobeam=True, add_foregrounds=True, apply_window=True, max_cached=1):
        self.noise_gen  = noise.NoiseGen(version=version,model=model,extract_region=extract_region,extract_region_shape=extract_region_shape,extract_region_wcs=extract_region_wcs,ncache=max_cached,verbose=False)
        self.signal_gen = signal.SignalGen(cmb_type=cmb_type, dobeam=dobeam, add_foregrounds=add_foregrounds, apply_window=apply_window, max_cached=max_cached, model=model, extract_region=extract_region,extract_region_shape=extract_region_shape, extract_region_wcs=extract_region_wcs)

    ### create wrapper for functions
    # get_fg
    # get_signal
    # get_noise
    # get_cmb
