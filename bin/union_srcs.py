from __future__ import print_function
from orphics import maps,io,cosmology,stats
from pixell import enmap,reproject,utils
import numpy as np
import os,sys
from soapack import interfaces as sints

"""
A brute-force way of finding duplicates among two catalogs.
This scales poorly with size of the catalogs.
"""

a = np.deg2rad([ [0,0 ],
               [1,1]])

print(np.rad2deg(utils.angdist(a,a[::-1]) ))
sys.exit()

def get_union(eras,edecs,bras,bdecs,dcut):

    # We first trim the big catalog
    dupes = 0
    for bra1,bdec1 in zip(bras,bdecs):
        ndupes = 0
        for bra2,bdec2 in zip(bras,bdecs):
            dra = bra1-bra2
            ddec = bdec1-bdec2
            dist = np.sqrt(dra**2.+ddec**2.)
            if dist*60.<dcut:
                dupes = dupes + 1
                ndupes = ndupes + 1
        if ndupes==0:
            fras.append(era)
            fdecs.append(edec)

    
    dupes = 0
    fras = []
    fdecs = []
    for era,edec in zip(eras,edecs):
        ndupes = 0
        for bra,bdec in zip(bras,bdecs):
            dra = era-bra
            ddec = edec-bdec
            dist = np.sqrt(dra**2.+ddec**2.)
            if dist*60.<dcut:
                dupes = dupes + 1
                ndupes = ndupes + 1
        print(ndupes)
        assert ndupes==0 or ndupes==1
        if ndupes==0:
            fras.append(era)
            fdecs.append(edec)
    fras = np.append(np.array(fras) , bras)
    fdecs = np.append(np.array(fdecs), bdecs)
                
    print(dcut,dupes)
    return fras,fdecs


eras,edecs,snrs,sizes = sints.get_act_mr3f_extended_sources(sn_cut=90.,size_cut=0.0)
print("Number of bright extended sources : ", len(eras))
bras,bdecs = sints.get_act_mr3f_cut_sources()
print("Number of bright cut sources : ", len(bras))

for dcut in [1.,0.5,0.1,0.05]:
    oras,odecs = get_union(eras,edecs,bras,bdecs,dcut)
    print(len(oras))

