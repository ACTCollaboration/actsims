import numpy as np
import os,sys
from pixell import enmap

map_root = "/home/msyriac/data/act/maps/mr3/"
mask_root = "/home/msyriac/data/act/maps/steve/"

class DataModel(object):
    def __init__(self,season,array,patch):
        self.mask_a = enmap.read_map("%s%s_mask_run_180323_master_apo_w0.fits" % (mask_root,patch))
        self.shape,self.wcs = self.mask_a.shape,self.mask_a.wcs
        self.region = patch
        self.array = array
        self.season = season
        self.freqs = {'pa1':['f150'],'pa2':['f150'],'pa3':['f090','f150']}[array]
        self.nfreqs = len(self.freqs)
        
    def get_inv_var(self):
        orets = []
        for freq in self.freqs:
            rets = []
            for k in range(4):
                pref = '_'.join([self.season,self.region,self.array,freq])
                rets.append( enmap.read_map("%s%s_nohwp_night_3pass_4way_set%d_ivar.fits" % (map_root,pref,k))[None] )
                # DEBUGGING rets.append( enmap.read_map("%s%s_nohwp_night_3pass_4way_set0_ivar.fits" % (map_root,pref))[None]*0+1. )

            orets.append(np.stack(rets))
        return enmap.enmap(np.stack(orets),self.wcs)
    
    def get_map(self):
        orets = []
        for freq in self.freqs:
            rets = []
            for k in range(4):
                pref = '_'.join([self.season,self.region,self.array,freq])
                rets.append( enmap.read_map("%s%s_nohwp_night_3pass_4way_set%d_map_srcfree.fits" % (map_root,pref,k)) )
            orets.append(np.stack(rets))
        return enmap.enmap(np.stack(orets),self.wcs)

class bin2D(object):
    def __init__(self, modrmap, bin_edges):
        self.centers = (bin_edges[1:]+bin_edges[:-1])/2.
        self.digitized = np.digitize(np.ndarray.flatten(modrmap), bin_edges,right=True)
        self.bin_edges = bin_edges
    def bin(self,data2d,weights=None):
        
        if weights is None:
            res = np.bincount(self.digitized,(data2d).reshape(-1))[1:-1]/np.bincount(self.digitized)[1:-1]
        else:
            res = np.bincount(self.digitized,(data2d*weights).reshape(-1))[1:-1]/np.bincount(self.digitized,weights.reshape(-1))[1:-1]
        return self.centers,res
    
def binned_power(pmap,modlmap,bin_edges):
    s = bin2D(modlmap,bin_edges)
    return s.bin(pmap)




