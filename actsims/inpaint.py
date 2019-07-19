from __future__ import print_function
from orphics import pixcov,io,cosmology,maps
from pixell import enmap,reproject
import numpy as np
from soapack import interfaces as sints
import os,sys
from enlib import bench

def inpaint_map_white(imap,ivar,fn_beam,union_sources_version=None,noise_pix = 20,hole_radius = 3.,plots=False):
    """

    Inpaints a map under the assumption of inhomogenous but white uncorrelated instrument noise.
    Pros: no products needed other than map-maker outputs of map and ivar, good inpainting despite crappy
    noise model.
    Cons: noise model has to be built for each source, so this is actually quite slow.

    imap -- (npol,Ny,Nx)
    ivar -- (Ny,Nx)
    fn_beam -- lambda ells: beam(ells)
    """

    ras,decs = sints.get_act_mr3f_union_sources(version=union_sources_version)
    cmb_theory_fn = lambda s,l: cosmology.default_theory().lCl(s,l)

    gtags = []
    gdicts = {}
    pcoords = []
    for i,(ra,dec) in enumerate(zip(ras,decs)):
        sel = reproject.cutout(ivar, ra=np.deg2rad(ra), dec=np.deg2rad(dec), pad=1, corner=False,npix=noise_pix,return_slice=True)
        if sel is None: continue
        civar = ivar[sel]
        if np.any(civar<=0): continue
        modrmap = civar.modrmap()
        modlmap = civar.modlmap()
        res = maps.resolution(civar.shape,civar.wcs)
        cimap = imap[sel]
        print("Inpainting: built noise model for source ",i," / ",len(ras))
        if plots: 
            for p in range(3): io.plot_img(cimap[p],os.environ['WORK']+"/cimap_%d_%s" % (p,str(i).zfill(2)))
            mimap = cimap.copy()
            mimap[...,modrmap<np.deg2rad(hole_radius/60.)] = np.nan
            for p in range(3): io.plot_img(mimap[p],os.environ['WORK']+"/masked_cimap_%d_%s" % (p,str(i).zfill(2)))
        
        scov = pixcov.scov_from_theory(modlmap,cmb_theory_fn,fn_beam,iau=False)
        ncov = pixcov.ncov_from_ivar(civar)
        pcov = scov + ncov
        gdicts[i] = pixcov.make_geometry(hole_radius=np.deg2rad(hole_radius/60.),n=noise_pix,deproject=True,iau=False,pcov=pcov,res=res)
        pcoords.append(np.array((dec,ra)))
        gtags.append(i)

    if len(gtags)>0: 
        pcoords = np.stack(pcoords).swapaxes(0,1)
        result = pixcov.inpaint(imap,pcoords,deproject=True,iau=False,geometry_tags=gtags,geometry_dicts=gdicts,verbose=True)

    if plots:
        for i,(ra,dec) in enumerate(zip(ras,decs)):
            sel = reproject.cutout(ivar, ra=np.deg2rad(ra), dec=np.deg2rad(dec), pad=1, corner=False,npix=noise_pix,return_slice=True)
            if sel is None: continue
            civar = ivar[sel]
            if np.any(civar<=0): continue
            modrmap = civar.modrmap()
            modlmap = civar.modlmap()
            res = maps.resolution(civar.shape,civar.wcs)
            cimap = result[sel]
            print("Inpainted ", ra,dec)
            if plots: 
                for p in range(3): io.plot_img(cimap[p],os.environ['WORK']+"/inpainted_cimap_%d_%s" % (p,str(i).zfill(2)))

    return result


def inpaint_map_const_cov(imap,mask,union_sources_version=None,noise_pix = 20,hole_radius = 3.,plots=False):
    """

    Inpaints a map under the assumption of constant 2D Fourier covariance. This uses the average PS
    of the full map for the noise model at each source location and thus does not handle inhomogenity.
    Pros: no products needed other than map-maker outputs of map
    Cons: 

    imap -- (npol,Ny,Nx)
    ivar -- (Ny,Nx)
    """

    ras,decs = sints.get_act_mr3f_union_sources(version=union_sources_version)
    kmap = enmap.fft(mask*imap,normalize='phys')

    gtags = []
    gdicts = {}
    pcoords = []
    for i,(ra,dec) in enumerate(zip(ras,decs)):
        sel = reproject.cutout(ivar, ra=np.deg2rad(ra), dec=np.deg2rad(dec), pad=1, corner=False,npix=noise_pix,return_slice=True)
        if sel is None: continue
        civar = ivar[sel]
        if np.any(civar<=0): continue
        modrmap = civar.modrmap()
        modlmap = civar.modlmap()
        res = maps.resolution(civar.shape,civar.wcs)
        cimap = imap[sel]
        print(ra,dec)
        if plots: 
            for p in range(3): io.plot_img(cimap[p],os.environ['WORK']+"/cimap_%d_%s" % (p,str(i).zfill(2)))
            mimap = cimap.copy()
            mimap[...,modrmap<np.deg2rad(hole_radius/60.)] = np.nan
            for p in range(3): io.plot_img(mimap[p],os.environ['WORK']+"/masked_cimap_%d_%s" % (p,str(i).zfill(2)))
        
        scov = pixcov.scov_from_theory(modlmap,cmb_theory_fn,fn_beam,iau=False)
        ncov = pixcov.ncov_from_ivar(civar)
        pcov = scov + ncov
        gdicts[i] = pixcov.make_geometry(hole_radius=np.deg2rad(hole_radius/60.),n=noise_pix,deproject=True,iau=False,pcov=pcov,res=res)
        pcoords.append(np.array((dec,ra)))
        gtags.append(i)

    if len(gtags)>0: 
        pcoords = np.stack(pcoords).swapaxes(0,1)
        result = pixcov.inpaint(imap,pcoords,deproject=True,iau=False,geometry_tags=gtags,geometry_dicts=gdicts,verbose=True)

    if plots:
        for i,(ra,dec) in enumerate(zip(ras,decs)):
            sel = reproject.cutout(ivar, ra=np.deg2rad(ra), dec=np.deg2rad(dec), pad=1, corner=False,npix=noise_pix,return_slice=True)
            if sel is None: continue
            civar = ivar[sel]
            if np.any(civar<=0): continue
            modrmap = civar.modrmap()
            modlmap = civar.modlmap()
            res = maps.resolution(civar.shape,civar.wcs)
            cimap = result[sel]
            print("Inpainted ", ra,dec)
            if plots: 
                for p in range(3): io.plot_img(cimap[p],os.environ['WORK']+"/inpainted_cimap_%d_%s" % (p,str(i).zfill(2)))

    return result
        

def test():
    # season = "s15"
    # patch = "deep56"
    # array = "pa3_f090"

    season = "s15"
    patch = "boss"
    array = "pa3_f150"


    dm = sints.ACTmr3()
    imap = dm.get_split(season,patch,array,0,ncomp=None,srcfree=True)
    ivar = dm.get_split_ivar(season,patch,array,0)
    fbeam = lambda x: dm.get_beam(x,season,patch,array)
    with bench.show("inpaint with plots"):
        inpainted = inpaint_map_white(imap,ivar,fbeam,plots=True)
    with bench.show("inpaint"):
        inpainted = inpaint_map_white(imap,ivar,fbeam,plots=False)



