import numpy as np, os
from orphics import io
import warnings

actsim_root = os.path.dirname(os.path.realpath(__file__))

try: paths = io.config_from_yaml(os.path.join(actsim_root, "../inputParams/paths.yml"))
except:
    paths = io.config_from_yaml(os.path.join(actsim_root, "../inputParams/paths_example.yml"))
    warnings.warn("No input/paths.yml found. Using version controlled input/paths_example.yml. Please copy and edit your local version.")

beam_root = paths['beam_root']

def __interpolate_beam__(l,f,l_interp):
    from scipy.interpolate import interp1d
    f_interp = interp1d(l, f, bounds_error=False, fill_value=(0,0))(l_interp)
    return (l_interp, f_interp)
    
def load_unnormalized_beam(season, array, freq, patch, l_interp=None):
    sapf = '{}_{}_{}_{}'.format(season, array, patch, freq)
    beam_file_temp = os.path.join(beam_root, 'mr3c_{}_{}_{}_nohwp_night_beam_tform_jitter_{}_181220.txt')
    beam_file      = beam_file_temp.format(season, array, freq, patch) 
    print("loading %s" %beam_file)
    #assert(os.path.exists(beam_file)) 
    l, f      = np.loadtxt(beam_file, unpack=True, usecols=[0,1])
    if l_interp is not None:
        l, f  = __interpolate_beam__(l,f,l_interp)
    return (l,f)

def load_normalized_beam(season, array, patch, freq, l_interp=None):
    l, f = load_unnormalized_beam(season, array, patch, freq, l_interp=None)
    assert(np.isclose(l[0],0))
    f    = f/f[0]
    if l_interp is not None:
        l, f  = __interpolate_beam__(l,f,l_interp)
    return (l,f)
    
def load_effective_beam(season, array, patch, freq, l_interp=None):
    raise NotImplemented()
