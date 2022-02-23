from __future__ import print_function
from orphics import pixcov,io,cosmology,maps
from pixell import enmap,reproject
import numpy as np
from soapack import interfaces as sints
import os,sys
from enlib import bench

def _get_temp_root_dir():
    try:
        return sints.dconfig['actsims']['temp_path']
    except KeyError:
        print("Error: Key temp_path not found in actsims section of ~/.soapack.yml. Please add this key and point it to the directory where you wish to hold temporary actsims files.")
        raise KeyError

def _get_radec_filename(rootdir):
    return rootdir + "/all_ra_dec.txt"

def _get_pcoords_filename(rootdir):
    return rootdir + "/pcoords.npy"

def _get_gtags_filename(rootdir):
    return rootdir + "/gtags.npy"

def _get_gdicts_filename(rootdir,key,item):
    return rootdir + "/gdicts_%d_%s.npy" % (key,item)

def load_cached_inpaint_geometries(cache_name):
    rootdir = _get_temp_root_dir() + cache_name
    assert os.path.exists(rootdir)
    
    ras,decs = np.loadtxt(_get_radec_filename(rootdir),unpack=True)
    pcoords = np.load(_get_pcoords_filename(rootdir))
    gtags = np.load(_get_gtags_filename(rootdir))
    gdicts = {}
    for key in gtags:
        gdicts[key] = {}
        for item in ['covsqrt','hole_radius','m1','m2','meanmul','ncomp','n','res']:
            gdicts[key][item] = np.load(_get_gdicts_filename(rootdir,key,item))


    return ras,decs,gtags,pcoords,gdicts

def save_cached_inpaint_geometries(cache_name,ras,decs,gtags,pcoords,gdicts):
    rootdir = _get_temp_root_dir() + cache_name
    assert not(os.path.exists(rootdir))
    os.mkdir(rootdir)
    
    io.save_cols(_get_radec_filename(rootdir),(ras,decs))
    np.save(_get_pcoords_filename(rootdir),pcoords)
    np.save(_get_gtags_filename(rootdir),gtags)

    for key in gdicts.keys():
        gd = gdicts[key]
        for item in gd.keys():
            np.save(_get_gdicts_filename(rootdir,key,item),gdicts[key][item])
    
    
def inpaint_map_white(imap,ivar,fn_beam,union_sources_version=None,noise_pix = 40,hole_radius = 6.,plots=False,cache_name=None,verbose=True,include_cmb=True,max_srcs=None,ppath=None):
    """

    Inpaints a map under the assumption of inhomogenous but white uncorrelated instrument noise.
    Pros: no products needed other than map-maker outputs of map and ivar, good inpainting despite crappy
    noise model.
    Cons: noise model has to be built for each source, so this is actually quite slow.

    imap -- (npol,Ny,Nx)
    ivar -- (Ny,Nx)
    fn_beam -- lambda ells: beam(ells)
    cache_name -- a unique string identifying the catalog+map/array/frequency/split combination to/from which the geometries are cached
    """
    if imap.ndim!=3:
        raise Exception("imap has to be of shape (ncomp,Ny,Nx)")
    ncomp = imap.shape[0]
   
    if cache_name is not None:
        cache_name = cache_name + "_catversion_%s" % union_sources_version
        try:
            ras,decs,gtags,pcoords,gdicts = load_cached_inpaint_geometries(cache_name)
            do_geoms = False
            if verbose: print("actsims.inpaint: loaded cached geometries for ", cache_name)
        except:
            if verbose: print("actsims.inpaint: no cached geometries found for ", cache_name, ". Generating and saving...")
            do_geoms = True
    else:
        do_geoms = True

    if do_geoms:
        iras,idecs = sints.get_act_mr3f_union_sources(version=union_sources_version)
        iras = iras[:max_srcs]
        idecs = idecs[:max_srcs]
        cmb_theory_fn = lambda s,l: cosmology.default_theory().lCl(s,l)
        gtags = []
        gdicts = {}
        pcoords = []

        ras = []
        decs = []
        for i,(ra,dec) in enumerate(zip(iras,idecs)):
            sel = reproject.cutout(ivar, ra=np.deg2rad(ra), dec=np.deg2rad(dec), pad=1, corner=False,npix=noise_pix,return_slice=True)
            if sel is None: continue
            ras.append(ra)
            decs.append(dec)

        for i,(ra,dec) in enumerate(zip(ras,decs)):
            sel = reproject.cutout(ivar, ra=np.deg2rad(ra), dec=np.deg2rad(dec), pad=1, corner=False,npix=noise_pix,return_slice=True)
            if sel is None: continue
            civar = ivar[sel]
            if np.any(civar<=0): continue
            modrmap = civar.modrmap()
            modlmap = civar.modlmap()
            res = maps.resolution(civar.shape,civar.wcs)
            cimap = imap[sel]
            if verbose: print("actsims.inpaint: building noise model for source ",i+1," / ",len(ras))
            if include_cmb:
                scov = pixcov.scov_from_theory(modlmap,cmb_theory_fn,fn_beam,iau=False,ncomp=ncomp)
            else:
                scov = 0
            ncov = pixcov.ncov_from_ivar(civar,ncomp=ncomp)
            pcov = scov + ncov
            gdicts[i] = pixcov.make_geometry(hole_radius=np.deg2rad(hole_radius/60.),n=noise_pix,deproject=True,iau=False,pcov=pcov,res=res)
            pcoords.append(np.array((dec,ra)))
            gtags.append(i)
        if len(gtags)>0: 
            pcoords = np.stack(pcoords).swapaxes(0,1)
        if cache_name is not None:
            save_cached_inpaint_geometries(cache_name,ras,decs,gtags,pcoords,gdicts)
            if verbose: print("actsims.inpaint: cached geometries for ",cache_name)

    from enlib import bench
    if len(gtags)>0: 
        with bench.show("inpaint"):
            result = pixcov.inpaint(imap,pcoords,deproject=True,iau=False,geometry_tags=gtags,geometry_dicts=gdicts,verbose=verbose)
    else:
        print("Nothing to inpaint")

    if ppath is None:
        ppath = os.environ['WORK']
    if plots:
        for i,(ra,dec) in enumerate(zip(ras,decs)):
            sel = reproject.cutout(ivar, ra=np.deg2rad(ra), dec=np.deg2rad(dec), pad=1, corner=False,npix=noise_pix,return_slice=True)
            if sel is None: continue
            civar = ivar[sel]
            if np.any(civar<=0): continue
            modrmap = civar.modrmap()
            modlmap = civar.modlmap()
            res = maps.resolution(civar.shape,civar.wcs)
            cimap = imap[sel]
            for p in range(ncomp): io.plot_img(cimap[p],ppath+"/cimap_%d_%s" % (p,str(i).zfill(2)))
            mimap = cimap.copy()
            mimap[...,modrmap<np.deg2rad(hole_radius/60.)] = np.nan
            for p in range(ncomp): io.plot_img(mimap[p],ppath+"/masked_cimap_%d_%s" % (p,str(i).zfill(2)))
            cimap = result[sel]
            for p in range(ncomp): io.plot_img(cimap[p],ppath+"/inpainted_cimap_%d_%s" % (p,str(i).zfill(2)))
            
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



