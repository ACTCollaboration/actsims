import numpy as np, os
from . import flipperDict

input_param_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../inputParams')
input_path      = lambda x : os.path.join(input_param_dir, x)

signal_dict     = flipperDict.flipperDict()
signal_dict.read_from_file(input_path('signal.dict'))

def __interpolate_beam__(l,f,l_interp):
    from scipy.interpolate import interp1d
    f_interp = interp1d(l, f, bounds_error=False, fill_value=(0,0))(l_interp)
    return (l_interp, f_interp)
    
def load_unnormalized_beam(patch, season, array, freq, l_interp=None):
    psaf = '{}_{}_{}_{}'.format(patch, season, array, freq)
    try:     
        assert(signal_dict['beamNames'].has_key(psaf))
    except:
        print("allowed inputs are %s" %str(signal_dict['beamNames']))
   
    beam_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),signal_dict['beamNames'][psaf])
    print("loading %s" %beam_file)
    l, f      = np.loadtxt(beam_file, unpack=True, usecols=[0,1])
    if l_interp is not None:
        l, f  = __interpolate_beam__(l,f,l_interp)
    return (l,f)

def load_normalized_beam(patch, season, array, freq, l_interp=None):
    l, f = load_unnormalized_beam(patch, season, array, freq, l_interp=None)
    assert(np.isclose(l[0],0))
    f    = f/f[0]
    if l_interp is not None:
        l, f  = __interpolate_beam__(l,f,l_interp)
    return (l,f)
    
def load_effective_beam(patch, season, array, freq, l_interp=None):
    raise NotImplemented()
