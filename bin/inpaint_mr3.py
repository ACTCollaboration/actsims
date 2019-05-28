from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
from orphics import maps,io,cosmology,pixcov
from pixell import enmap,reproject
import numpy as np
import os,sys
import soapack.interfaces as sints
from enlib import bench

apatch = sys.argv[1]

out_dir = "/scratch/r/rbond/msyriac/data/depot/actsims/inpainted/"

plot_img = lambda x,y: io.plot_img(x,os.environ['WORK']+"/new_mr3f/"+y)
noise_pix = 60
cmb_theory_fn = lambda s,l: cosmology.default_theory().lCl(s,l)
hole_radius = 9.
res = np.deg2rad(0.5/60.)

def plot_cutout(nsplits,cutout,ivars,tag="",skip_plots=False):
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
            if not(skip_plots): plot_img(cutout[s,p],"cutout_%s.png" % sname)
            masked = cutout[s,p].copy()
            modrmap = masked.modrmap()
            cutradius = np.deg2rad(hole_radius/60.)
            masked[modrmap<cutradius] = np.nan
            if not(skip_plots): plot_img(masked,"masked_%s.png" % sname)
        masked = ivars[s,0].copy()
        masked[modrmap<cutradius] = np.nan
        iname = "%s%s_%s_%s_split_%d_%d" % (tag,season,patch,array,s,sid)
        if not(skip_plots): plot_img(masked,"ivars_masked_%s.png" % iname)
        retvars[s,0,modrmap<cutradius] = retvars[s,0,modrmap>=cutradius].mean()
        if not(skip_plots): plot_img(retvars[s,0],"ivars_filled_%s.png" % iname)
    return retvars

ras,decs = np.loadtxt("/home/r/rbond/sigurdkn/project/actpol/maps/mr3f_20190502/cat_bright_tot.txt",unpack=True,usecols=[0,1])

dm = sints.ACTmr3(calibrated=False)

for season in dm.seasons:
    for patch in dm.patches:
        if patch!=apatch: continue
        for array in dm.arrays:
            if array!="pa3_f090": continue # !!!
            try: 
                splits = dm.get_splits(season,patch,array,ncomp=None,srcfree=True)[0]
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
                if (cutout is not None): 
                    sel_ivar = reproject.cutout(ivars, ra=np.deg2rad(ra), dec=np.deg2rad(dec), pad=1, corner=False,npix=noise_pix,return_slice=True)
                    cut_ivar = ivars[sel_ivar]
                    retv = plot_cutout(nsplits,cutout,cut_ivar,skip_plots=True)
                    ivars[sel_ivar] = retv
                    ids.append(sid)
                    sel_ivars.append(sel_ivar)


            # Save ivar map
            #....
            for i in range(nsplits):
                fname = out_dir+os.path.basename(dm.get_split_ivar_fname(season,patch,array,i))
                enmap.write_map(fname,ivars[i])
            
            # Inpaint each split
            beam_fn = lambda x: dm.get_beam(season=season,patch=patch,array=array,ells=x)
            inpainted = []
            for i in range(nsplits):
                gtags = []
                gdicts = {}
                pcoords = []
                for sindex,sid in enumerate(ids):
                    cutout = reproject.cutout(splits[i], ra=np.deg2rad(ras[sid]), dec=np.deg2rad(decs[sid]), pad=1, corner=False,npix=noise_pix)
                    if np.std(cutout)<1e-3: continue
                    with bench.show("geometry"):
                        pcov = pixcov.pcov_from_ivar(noise_pix,np.deg2rad(decs[sid]),np.deg2rad(ras[sid]),ivars[i,0],cmb_theory_fn,beam_fn,iau=False)
                        gdicts[sid]  = pixcov.make_geometry(hole_radius=np.deg2rad(hole_radius/60.),n=noise_pix,deproject=True,iau=False,pcov=pcov,res=res)
                        gtags.append(sid)
                    pcoords.append(np.array((decs[sid],ras[sid])))
                if len(gtags)>0: 
                    pcoords = np.stack(pcoords).swapaxes(0,1)
                    result = pixcov.inpaint(splits[i],pcoords,deproject=True,iau=False,geometry_tags=gtags,geometry_dicts=gdicts,verbose=True)
                else:
                    result = splits[i].copy()
                inpainted.append( result)
            inpainted = enmap.enmap(np.stack(inpainted),wcs)

                
            # Verify inpainting
            for sid,(ra,dec) in enumerate(zip(ras,decs)):
                cutout = reproject.cutout(inpainted, ra=np.deg2rad(ra), dec=np.deg2rad(dec), pad=1, corner=False,npix=noise_pix)
                if (cutout is not None): 
                    sel_ivar = reproject.cutout(ivars, ra=np.deg2rad(ra), dec=np.deg2rad(dec), pad=1, corner=False,npix=noise_pix,return_slice=True)
                    cut_ivar = ivars[sel_ivar]
                    retv = plot_cutout(nsplits,cutout,cut_ivar,tag="inpainted_")


            # Save inpainted map
            #....

            for i in range(nsplits):
                fname = out_dir+os.path.basename(dm.get_split_fname(season,patch,array,i,srcfree=True))
                enmap.write_map(fname,inpainted[i])
                
