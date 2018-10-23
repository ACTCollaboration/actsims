
import scipy, enlib.enmap, numpy as np

def mySmoothing(inmap, N = 100):

    smoothed = scipy.ndimage.filters.gaussian_filter(inmap, N)
    
    new = np.zeros(smoothed.shape)
    new[smoothed > .9] = 1.

    smoothedagain = scipy.ndimage.filters.gaussian_filter(new, N)

    return enlib.enmap.ndmap(smoothedagain, inmap.wcs)
