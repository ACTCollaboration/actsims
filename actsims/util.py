#-
# util.py
#-
import os, numpy as np
from collections import OrderedDict

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
        has_it = nested_dict.has_key(keys[0])
        return has_key(nested_dict[keys[0]], keys[1:]) if has_it else False
    else:
        return nested_dict.has_key(keys[0])

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
