from __future__ import print_function
from orphics import maps,io,cosmology,stats
from pixell import enmap
import numpy as np
import os,sys
from szar import foregrounds as fg

nu = 150.
ells = np.arange(0,8000,1)
theory = cosmology.default_theory()
cltt = theory.lCl('TT',ells)
clyy = fg.power_tsz(ells,nu,nu2=None,A_tsz=None,fill_type="extrapolate",tcmb=None)

totres = fg.power_cibp(ells,nu,nu2=None) + fg.power_cibc(ells,nu,nu2=None) + fg.power_radps(ells,nu,nu2=None) + fg.power_ksz_reion(ells) + + fg.power_ksz_late(ells)

def resmodel(ells,a1,a2,e1,e2,ell0=3000):
    plaw = lambda a,e: a * (ells/ell0)**(e)
    return plaw(a1,e1) + plaw(a2,e2)

clres = resmodel(ells,1e-5,1e-5,0.5,1.0)

pl = io.Plotter(xyscale='linlog',scalefn=lambda x: x**2./2./np.pi)
pl.add(ells,cltt,lw=3,color='k')
pl.add(ells,clyy)
#pl.add(ells,clres)
pl.add(ells,totres,ls='--')
pl.done()


"""

For each freq1,freq2:
res = ps - cmb - tsz(freq1,freq2)
Fit a1*ell**e1 + a2*ell**e2 to res
With 2000 < ell < 8000 for act x act
With 2000 < ell < 6000 for act x hfi
With 2000 < ell < 3000 for act x lfi
With 2000 < ell < 6000 for hfi x hfi
With 2000 < ell < 3000 for lfi x lfi

cov(freq1,freq2) 

"""
