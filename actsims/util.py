#-
# util.py
#-
import os,sys, numpy as np
from collections import OrderedDict
from soapack import interfaces as dmint

class _SeedTracker(object):
    def __init__(self):
        self.CMB     = 0
        self.FG      = 1
        self.PHI     = 2
        self.NOISE   = 3
        self.POISSON = 4
        self.COMPTONY = 5

        #quick-srcfree is maxmally correlated with 15mJy sims
        self.fgdict  = {'15mjy': 0, '100mjy': 1, 'srcfree': 2, 'quick-srcfree':0,'comptony': 3}

        self.dmdict  = {'act_mr3':0,'act_c7v5':1,'planck_hybrid':2}

    def get_cmb_seed(self, set_idx, sim_idx):
        return (set_idx, 0, self.CMB, sim_idx)

    def get_fg_seed(self, set_idx, sim_idx, fg_type):
        assert(fg_type in self.fgdict.keys())
        return (set_idx, 0, self.FG, sim_idx, self.fgdict[fg_type])

    def get_phi_seed(self, set_idx, sim_idx):
        return (set_idx, 0, self.PHI, sim_idx)

    def get_noise_seed(self, set_idx, sim_idx, data_model, season, patch, array, patch_id=None):
        ret = (set_idx, 0, self.NOISE, sim_idx)
        dm  = data_model
        
        assert(dm.name in self.dmdict.keys())
        sid =  dm.seasons.index(season) if dm.seasons is not None else 0
        pid =  dm.patches.index(patch) if dm.patches is not None else 0
        aid = list(dm.array_freqs.keys()).index(array)
        ret = ret + (self.dmdict[dm.name],sid,pid,aid)
        if patch_id is not None:
            ret = ret + (patch_id,)
        return ret

    def get_poisson_seed(self, set_idx, sim_idx):
        return (set_idx, 0, self.POISSON, sim_idx)

seed_tracker = _SeedTracker()


#########################


def is_empty(s):
    '''
        if s is either '' or None, this function returns True. Otherwise, false.
    '''
    return True if (s is None) or (len(s) == 0) else False

def create_dict(*idxes):
    '''
        create nested dictionary with the given idxes
    '''

    height  = len(idxes)
    output = {}
    
    stack  = []
    stack.append(output)

    for depth in range(height):
        stack_temp = []
        while len(stack) > 0:
            cur_elmt = stack.pop()
            for idx in idxes[depth]:
                cur_elmt[idx] = {}
                stack_temp.append(cur_elmt[idx])
        stack = stack_temp

    return output

def get_from_dict(nested_dict, keys, safe=True):
    if not type(keys) is tuple: keys = (keys,)
    if safe and not has_key(nested_dict, keys): return None

    if(len(keys) > 1):
        return get_from_nested_dict(nested_dict[keys[0]], keys[1:], False)
    else:
        return nested_dict[keys[0]]


def has_key(nested_dict, keys):
    ''' search through nested dictionary to fine the elements '''
    if not type(keys) is tuple: keys = (keys,)
    if not type(nested_dict) == dict: return False

    if(len(keys) > 1):
        has_it = keys[0] in nested_dict
        return has_key(nested_dict[keys[0]], keys[1:]) if has_it else False
    else:
        return keys[0] in nested_dict

def get_num_leaves(nested_dict):
    ''' search how many elements are in nested dictionary '''
    count = 0
    for key in nested_dict.keys():
        if isinstance(nested_dict[key], dict):
            count += get_num_elements(nested_dict[key])
        else:
            count += 1
    return count


# decorators
def static_vars(**kwargs):
    '''https://stackoverflow.com/questions/279561/what-is-the-python-equivalent-of-static-variables-inside-a-function'''
    def decorate(func):
        for k in kwargs:
            setattr(func, k, kwargs[k])
        return func
    return decorate

def create_dir(path_to_dir):
    ''' check wether the directory already exists. if not, create it '''
    exit = 0       # exit code

    if not os.path.isdir(path_to_dir):
        os.makedirs(path_to_dir)
        exit = 0
    else:
        exit = 1

    return exit

# took it from quicklens/util.py
class memorize(object):
    """ a simple memoize decorator (http://en.wikipedia.org/wiki/Memoization) """
    def __init__(self, func):
        self.func = func
        self.cache = {}
    def __call__(self, *args):
        if args in self.cache:
            print("returning from cache.")
            return self.cache[args]
        else:
            v = self.func(*args)
            self.cache[args] = v
            return v

def mkdir(dirpath,comm=None):
    if comm is None:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
    exists = os.path.exists(dirpath)
    comm.Barrier()
    if comm.Get_rank()==0: 
        if not (exists):
            os.makedirs(dirpath)
    return exists
