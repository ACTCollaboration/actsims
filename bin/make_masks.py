from __future__ import print_function
from orphics import maps,io,cosmology
from pixell import enmap,fft
import numpy as np
import os,sys
from soapack import interfaces as sints

"""
This script produces analysis masks (binary mask + apodization)
for each region that will be used in noise template generation,
1-tile ILC and 1-tile lensing. It simply takes the apodized
binary masks from Steve, pads it a bit and then pads it
a bit more to make it FFT friendly.
"""


out_version = "padded_v1"
in_versions = {'actpol': "mr3c_20190215_pickupsub_190301",'advact': "mr3c_20190215_pickupsub_190303"}
patches = {'actpol':['deep8','deep1','deep5','deep6','deep56','boss'],'advact':['patch00%d' % i for i in range(9)]}
pads = {'deep':200, 'boss':400, 'patch': 600}

#out_version = "padded_deep_v2"
#in_versions = {'actpol': "180323"}
#patches = {'actpol':['deep1','deep5','deep6','deep56','boss']}
#pads = {'deep':1, 'boss':1}


out_path = "/home/r/rbond/msyriac/data/act/maps/steve/%s/" % out_version


for survey in in_versions.keys():
    for patch in patches[survey]:

        if 'deep' in patch:
            pad = pads['deep']
        elif patch=='boss':
            pad = pads['boss']
        elif 'patch' in patch:
            pad = pads['patch']
        else:
            raise ValueError

        mpatch = patch
        mask = sints.get_act_mr3_crosslinked_mask(mpatch,
                                                  version=in_versions[survey],
                                                  kind='binary_apod',
                                                  season="s16" if survey=='advact' else None,array=None,
                                                  pad=pad)
        print(survey,patch)
        # FFT friendliness
        Ny,Nx = mask.shape[-2:]
        dNy = fft.fft_len(Ny,"above")
        dNx = fft.fft_len(Nx,"above")
        pny = dNy - Ny
        pnx = dNx - Nx
        pady1 = pny//2
        pady2 = pny - pady1
        padx1 = pnx//2
        padx2 = pnx - padx1
        mask = enmap.pad(mask,((pady1,padx1),(pady2,padx2)))
        assert mask.shape[-2]==dNy
        assert mask.shape[-1]==dNx

        enmap.write_map(out_path+"%s.fits" % patch,mask)
        io.hplot(enmap.downgrade(mask,8),out_path+"%s" % patch)
        
