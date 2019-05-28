from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
from orphics import maps,io,cosmology,pixcov
from pixell import enmap,reproject
import numpy as np
import os,sys
import soapack.interfaces as sints
from enlib import bench
from actsims import noise

# apatch = sys.argv[1]
# season = sys.argv[2]
# array = sys.argv[3]

dm = sints.ACTmr3(calibrated=False)

for season in ['s13','s14','s15']:
    for apatch in ['deep1','deep5','deep6','deep8','deep56','boss']:
        for array in ['pa1_f150','pa2_f150','pa3_f150','pa3_f090']:

            fname = "%s_%s_%s" % (season,apatch,array)
            try:
                splits = dm.get_splits(season=season,patch=apatch,arrays=[array],ncomp=1,srcfree=False)[0,:,0,...]
            except:
                continue
            ivars = dm.get_splits_ivar(season=season,patch=apatch,arrays=[array],ncomp=None)[0,:,0,...]
            cmap,_ = noise.get_coadd(splits,ivars,axis=0)

            io.hplot(cmap,os.environ['WORK']+"/new_mr3f/coadd_%s.png" % fname,min=-300,max=600,grid=True)
            splits = dm.get_splits(season=season,patch=apatch,arrays=[array],ncomp=1,srcfree=True)[0,:,0,...]
            cmap,_ = noise.get_coadd(splits,ivars,axis=0)
            io.hplot(cmap,os.environ['WORK']+"/new_mr3f/coadd_srcfree_%s.png" % fname,min=-300,max=600,grid=True)
sys.exit()

out_dir = "/scratch/r/rbond/msyriac/data/depot/actsims/inpainted/"

plot_img = lambda x,y,**kwargs: io.plot_img(x,os.environ['WORK']+"/new_mr3f/"+y,cmap='gray',**kwargs)
noise_pix = 60
cmb_theory_fn = lambda s,l: cosmology.default_theory().lCl(s,l)
hole_radius = 9.
res = np.deg2rad(0.5/60.)

def plot_cutout(nsplits,cutout,pcutout,ivars,tag="",skip_plots=False):
    pols = ['I','Q','U']
    retvars = ivars.copy()
    for s in range(nsplits):
        if np.std(cutout[s])<1e-3: 
            print("Skipping split %d as it seems empty" % s)
            continue
        iname = "%s%s_%s_%s_split_%d_%d" % (tag,season,patch,array,s,sid)
        if not(skip_plots): plot_img(ivars[s,0],"ivars_%s.png" % iname)
        for p in range(3):
            pol = pols[p]
            sname = "%s%s_%s_%s_split_%d_%s_%d" % (tag,season,patch,array,s,pol,sid)
            if not(skip_plots): 
                plot_img(cutout[s,p],"cutout_%s.png" % sname,lim=1000)
                plot_img(pcutout[s,p],"pcutout_%s.png" % sname,lim=1000)
            masked = cutout[s,p].copy()
            modrmap = masked.modrmap()
            cutradius = np.deg2rad(hole_radius/60.)
            masked[modrmap<cutradius] = np.nan
            #if not(skip_plots): plot_img(masked,"masked_%s.png" % sname)
        masked = ivars[s,0].copy()
        masked[modrmap<cutradius] = np.nan
        iname = "%s%s_%s_%s_split_%d_%d" % (tag,season,patch,array,s,sid)
        #if not(skip_plots): plot_img(masked,"ivars_masked_%s.png" % iname)
        retvars[s,0,modrmap<cutradius] = retvars[s,0,modrmap>=cutradius].mean()
        #if not(skip_plots): plot_img(retvars[s,0],"ivars_filled_%s.png" % iname)
    return retvars

ras,decs = np.loadtxt("inputParams/mr3_cut_sources.txt",unpack=True)

dm = sints.ACTmr3(calibrated=False)

for season in dm.seasons:
    if season == "s16": continue
    for patch in dm.patches:
        if patch!=apatch: continue
        for array in dm.arrays:
            try: 
                splits = dm.get_splits(season,patch,array,ncomp=None,srcfree=True)[0]
                psplits = dm.get_splits(season,patch,array,ncomp=None,srcfree=False)[0]
                print("Found %s %s %s" % (season,patch,array))
            except: 
                continue
            ivars = dm.get_splits_ivar(season,patch,array,ncomp=None)[0]
            wcs = ivars.wcs
            nsplits = ivars.shape[0]

            ids = []
            sel_ivars = []
            for sid,(ra,dec) in enumerate(zip(ras,decs)):
                cutout = reproject.cutout(splits, ra=np.deg2rad(ra), dec=np.deg2rad(dec), pad=1, corner=False,npix=noise_pix)
                pcutout = reproject.cutout(psplits, ra=np.deg2rad(ra), dec=np.deg2rad(dec), pad=1, corner=False,npix=noise_pix)
                if (cutout is not None): 
                    sel_ivar = reproject.cutout(ivars, ra=np.deg2rad(ra), dec=np.deg2rad(dec), pad=1, corner=False,npix=noise_pix,return_slice=True)
                    cut_ivar = ivars[sel_ivar]
                    retv = plot_cutout(nsplits,cutout,pcutout,cut_ivar,skip_plots=False)
                    ivars[sel_ivar] = retv
                    ids.append(sid)
                    sel_ivars.append(sel_ivar)


