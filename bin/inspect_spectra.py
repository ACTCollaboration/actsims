from __future__ import print_function
from orphics import maps,io,cosmology
from pixell import enmap
import numpy as np
import os,sys

import argparse
# Parse command line
parser = argparse.ArgumentParser(description='Do a thing.')
parser.add_argument("region", type=str,help='boss or deep56?')
parser.add_argument("--fft", action='store_true',help='A flag.')
args = parser.parse_args()

fftstr = "_fft" if args.fft else ""


def load_spec(a1,a2):
    def _load_spec(a1,a2):
        try:
            aseason1,aarray1 = a1
            str1 = "%s_%s" % (aseason1,aarray1)
        except:
            pfreq1 = a1
            str1 = "%s" % pfreq1
        try:
            aseason2,aarray2 = a2
            str2 = "%s_%s" % (aseason2,aarray2)
        except:
            pfreq2 = a2
            str2 = "%s" % pfreq2
        fname = "/scratch/r/rbond/msyriac/data/depot/actsims/spectra/spec%s_%s_%s_%s.txt" % (fftstr,args.region,str1,str2)
        ls,Cls = np.loadtxt(fname,unpack=True)
        return ls,Cls
    try:
        return _load_spec(a1,a2)
    except:
        return _load_spec(a2,a1)

pfreqs = ['030','044','070','100','143','217','353','545']
lfis = ['030','044','070']
hfis = ['100','143','217','353','545']
#pfreqs = lfis
nfreqs = len(pfreqs)

pl = io.Plotter(xyscale='linlog',scalefn = lambda x: x**2./np.pi,xlabel='$\\ell$',ylabel='$D_{\\ell}$')
for i in range(nfreqs):
    for j in range(i,nfreqs):
        ls,Cls = load_spec(pfreqs[i],pfreqs[j])

        if (pfreqs[i] in lfis)  or (pfreqs[j] in lfis):
            Cls = Cls[ls<600]
            ls = ls[ls<600]
        else:
            Cls = Cls[ls<3000]
            ls = ls[ls<3000]


        

        pl.add(ls,Cls,label = "%s x %s" % (pfreqs[i],pfreqs[j]))
pl._ax.set_ylim(1,1e7)
pl.hline(y=0)
pl.done(os.environ['WORK']+"/pspec%s.png" % fftstr)

