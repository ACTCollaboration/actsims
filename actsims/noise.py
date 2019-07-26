import numpy as np
import os,sys
from pixell import enmap,enplot,fft as pfft
from soapack import interfaces as sints
from actsims import util
from enlib import bench
import warnings
if 'fftw' not in pfft.engine: warnings.warn("No pyfftw found. Using much slower numpy fft engine.")

class NoiseGen(object):
    def __init__(self,version,model="act_mr3",extract_region=None,extract_region_shape=None,extract_region_wcs=None,ncache=1,verbose=False):
        """
        version: The version identifier for the filename of covsqrts on disk
        model: The name of an implemented soapack datamodel
        extract_region: An optional map whose footprint on to which the sims are made
        extract_region_shape: Instead of passing a map for extract_region, one can pass its shape and wcs
        extract_region_wcs: Instead of passing a map for extract_region, one can pass its shape and wcs
        ncache: The number of 

        """
        self._version = version
        self._model = model
        self.ncache = ncache
        self._ccache = {}
        self._icache = {}
        self.dm = sints.models[model](region=extract_region)
        self.verbose = verbose

    def save_covsqrt(self,covsqrt,season=None,patch=None,array=None,coadd=True,mask_patch=None):
        pout,cout,sout = get_save_paths(self._model,self._version,
                                        coadd=coadd,season=season,patch=patch,array=array,
                                        overwrite=False,mask_patch=mask_patch)
        fpath = "%s_covsqrt.fits" % (cout)
        enmap.write_map(fpath ,covsqrt)
        print(fpath,covsqrt.shape)

    def save_filter_noise(self,n2d,season=None,patch=None,array=None,coadd=True,mask_patch=None):
        pout,cout,sout = get_save_paths(self._model,self._version,
                                        coadd=coadd,season=season,patch=patch,array=array,
                                        overwrite=False,mask_patch=mask_patch)
        fpath = "%s_filter_noise.fits" % (cout)
        enmap.write_map(fpath ,n2d)
        print(fpath,n2d.shape,n2d.wcs)

    def load_covsqrt(self,season=None,patch=None,array=None,coadd=True,mask_patch=None,get_geometry=False):
        pout,cout,sout = get_save_paths(self._model,self._version,coadd=coadd,
                                        season=season,patch=patch,array=array,
                                        overwrite=False,mask_patch=mask_patch)
        fpath = "%s_covsqrt.fits" % (cout)
        if get_geometry: return enmap.read_map_geometry(fpath)
        ikey = '_'.join([str(x) for x in [season,patch,self.dm.array_freqs[array]]])
        try:
            covsqrt = self._ccache[fpath]
            ivars = self._icache[ikey]
            if self.verbose: print("Loaded cached covsqrt and ivars.")
        except:
            if self.verbose: print("Couldn't find covsqrt and ivars in cache. Reading from disk...")
            ivars = self.dm.get_splits_ivar(season=season,patch=patch,arrays=self.dm.array_freqs[array])
            covsqrt = enmap.read_map(fpath)
            print(fpath,covsqrt.shape)
            if len(self._ccache.keys())<self.ncache: 
                self._ccache[fpath] = covsqrt
                self._icache[ikey] = ivars
        return covsqrt,ivars
        
    def generate_sim(self,season=None,patch=None,array=None,seed=None,mask_patch=None,apply_ivar=True,ext_signal=0):
        covsqrt,ivars = self.load_covsqrt(season=season,patch=patch,array=array,mask_patch=mask_patch)
        sims,ivars = generate_noise_sim(covsqrt,ivars,seed=seed,dtype=self.dm.dtype)
        if apply_ivar: 
            sims = apply_ivar_window(ext_signal+sims,ivars)
            return sims
        else:
            return sims,ivars


    def save_sims(self,seed,sims,season,patch,array,mask_patch,coadd=True):
        pout,cout,sout = get_save_paths(self._model,self._version,coadd=coadd,
                                        season=season,patch=patch,array=array,
                                        overwrite=False,mask_patch=mask_patch)
        insplits = self.dm.get_nsplits(season,patch,array)
        freqs = self.dm.array_freqs[array]
        nfreqs,nsplits,npol,Ny,Nx = sims.shape
        assert nsplits == insplits
        assert len(freqs) == nfreqs
        for i in range(nfreqs):
            iarray = freqs[i]
            for j in range(nsplits):
                fname = sout+os.path.basename(self.dm.get_split_fname(season,patch,iarray,j,srcfree=True)).replace(".fits","_seed_%d.fits" % seed)
                enmap.write_map(fname,sims[i,j,:,:,:])
        

def apply_ivar_window(imaps,ivars):
    # input has shape (nfreqs,nsplits,npol,Ny,Nx)
    # Sanitize by thresholding and binary masking
    assert imaps.ndim==5
    nfreqs = imaps.shape[0]
    nsplits = imaps.shape[1]
    for ifreq in range(nfreqs):
        for isplit in range(nsplits):
            win = ivars[ifreq,isplit,0,...]
            imaps[ifreq,isplit,:,win==0.] = 0 # FIXME: make 0 comparison robust
    return imaps

def get_save_paths(model,version,coadd,season=None,patch=None,array=None,mkdir=False,overwrite=False,mask_patch=None):
    paths = sints.dconfig['actsims']

    try: assert paths['plot_path'] is not None
    except: paths['plot_path'] = "./"
    assert paths['covsqrt_path'] is not None
    assert paths['trial_sim_path'] is not None

    # Prepare output dirs
    pdir = "%s/%s/" % (paths['plot_path'] ,version) 
    cdir = "%s/%s/" % (paths['covsqrt_path'] ,version)
    sdir = "%s/%s/" % (paths['trial_sim_path'] ,version)
    
    if mkdir:
        exists1 = util.mkdir(pdir)
        exists2 = util.mkdir(cdir)
        exists3 = util.mkdir(sdir)
        if any([exists1,exists2,exists3]): 
            if not(overwrite): raise IOError
            warnings.warn("Version directory already exists. Overwriting.")

    if model=='planck_hybrid': 
        assert season is None
        suff = '_'.join([model,patch,array,"coadd_est_"+str(coadd)])
    else:
        suff = '_'.join([model,season,patch,array,"coadd_est_"+str(coadd)])


    pout = pdir + suff
    cout = cdir + suff
    sout = sdir

    if mask_patch is not None:
        if mask_patch != patch:
            pout = pout+"_"+mask_patch
            cout = cout+"_"+mask_patch
            sout = sout+mask_patch+"_"

    return pout,cout,sout


def get_n2d_data(splits,ivars,mask_a,coadd_estimator=False,flattened=False,plot_fname=None,dtype=None):
    assert np.all(np.isfinite(splits))
    assert np.all(np.isfinite(ivars))
    assert np.all(np.isfinite(mask_a))

    if coadd_estimator:
        coadd,civars = get_coadd(splits,ivars,axis=1)
        data  = splits - coadd[:,None,...]
        del coadd
    else:
        data = splits

    assert np.all(np.isfinite(data))
    if flattened:
        # sivars = np.sqrt(ivars)
        sivars   = ((1./ivars) - (1./civars[:,None,...]))**-0.5
        sivars[~np.isfinite(sivars)] = 0
        ffts = enmap.fft(data*mask_a*sivars,normalize="phys")
        if plot_fname is not None: plot(plot_fname+"_fft_maps",data*mask_a*sivars,quantile=0)
        wmaps = mask_a + enmap.zeros(ffts.shape,mask_a.wcs,dtype=dtype) # WARNING: type
        del ivars, data, splits
    else:
        assert np.all(np.isfinite(data*mask_a*ivars))
        ffts = enmap.fft(data*mask_a*ivars,normalize="phys")
        if plot_fname is not None: plot(plot_fname+"_fft_maps",data*mask_a*ivars,quantile=0)
        wmaps = ivars * mask_a
        del ivars, data, splits
    n2d = get_n2d(ffts,wmaps,coadd_estimator=coadd_estimator,plot_fname=plot_fname,dtype=dtype)
    assert np.all(np.isfinite(n2d))
    return n2d


def generate_noise_sim(covsqrt,ivars,seed=None,dtype=None):
    """
    Supports only two cases
    1) nfreqs>=1,npol=3
    2) nfreqs=1,npol=1
    """
    if isinstance(seed,int): seed = (seed,)
    assert np.all(np.isfinite(covsqrt))
    shape,wcs = covsqrt.shape,covsqrt.wcs
    Ny,Nx = shape[-2:]
    ncomps = covsqrt.shape[0]
    assert ncomps==covsqrt.shape[1]
    assert ((ncomps % 3) == 0) or (ncomps==1)
    nfreqs = 1 if ncomps==1 else ncomps // 3
    if ncomps==1: npol = 1
    else: npol = 3
    wmaps = enmap.extract(ivars,shape[-2:],wcs)
    nsplits = wmaps.shape[1]

    if dtype is np.float32: ctype = np.complex64 
    elif dtype is np.float64: ctype = np.complex128 

    # Old way with loop
    kmap = []
    for i in range(nsplits):
        if seed is None:
            np.random.seed(None)
        else:
            np.random.seed(seed+(i,))
        rmap = enmap.rand_gauss_harm((ncomps, Ny, Nx),covsqrt.wcs).astype(ctype)
        kmap.append( enmap.map_mul(covsqrt, rmap) )
    del covsqrt, rmap
    kmap = enmap.enmap(np.stack(kmap),wcs)
    outmaps = enmap.ifft(kmap, normalize="phys").real
    del kmap

    # Need to test this more ; it's only marginally faster and has different seed behaviour
    # covsqrt = icovsqrt 
    # np.random.seed(seed)
    # rmap = enmap.rand_gauss_harm((nsplits,ncomps,Ny, Nx),covsqrt.wcs)
    # kmap = enmap.samewcs(np.einsum("abyx,cbyx->cayx", covsqrt, rmap),rmap)
    # outmaps = enmap.ifft(kmap, normalize="phys").real

    # isivars = 1/np.sqrt(wmaps)
    isivars   = ((1./wmaps) - (1./wmaps.sum(axis=1)[:,None,...]))**0.5
    isivars[~np.isfinite(isivars)] = 0
    
    assert np.all(np.isfinite(outmaps))
    # Divide by hits
    for ifreq in range(nfreqs):
        outmaps[:,ifreq*npol:(ifreq+1)*npol,...] = outmaps[:,ifreq*npol:(ifreq+1)*npol,...] * isivars[ifreq,...] *np.sqrt(nsplits)

    retmaps = outmaps.reshape((nsplits,nfreqs,npol,Ny,Nx)).swapaxes(0,1)
    assert np.all(np.isfinite(retmaps))
    return retmaps,wmaps

    
def get_coadd(imaps,wts,axis):
    # sum(w*m)/sum(w)
    twt = np.sum(wts,axis=axis)
    retmap = np.nan_to_num(np.sum(wts*imaps,axis=axis)/twt)
    #retmap[np.isinf(np.abs(retmap))] = 0
    return retmap,twt



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

def null_off_diagonals(cov):
    ocov = cov.copy()
    for i in range(cov.shape[0]):
        for j in range(cov.shape[0]):
            if i==j: continue
            ocov[i,j] = 0
    return ocov
    
def null_pol_off_diagonals(cov):
    ocov = cov.copy()
    for i in range(cov.shape[0]):
        for j in range(cov.shape[0]):
            if i==j or (i==0 and j==3) or (i==3 and j==0): continue
            ocov[i,j] = 0
    return ocov
    
def get_covsqrt(ps,method="arrayops"):
    if method=="multipow":
        covsq = enmap.multi_pow(ps.copy(),0.5)
    elif method=="arrayops":
        from enlib import array_ops
        covsq = array_ops.eigpow(ps.copy(),0.5,axes=[0,1])
    covsq[:,:,ps.modlmap()<2] = 0
    assert np.all(np.isfinite(covsq))
    return enmap.enmap(covsq,ps.wcs)
    

def naive_power(f1,w1,f2=None,w2=None):
    if f2 is None:
        assert w2 is None
        f2 = f1
        w2 = w1
    norm = np.mean(w1*w2,axis=(-2,-1),keepdims=True)
    ret = np.real(f1*f2.conj()) / norm
    return ret
    
def noise_power(kmaps1,weights1,kmaps2=None,weights2=None,
                coadd_estimator=False,pfunc=None):
    # All weights must include the apodized mask
    # No weights must be applied before hand
    # coadds1 = (Ny,Nx) -- the coadd map
    # cweights1 = (Ny,Nx) -- the coadd weight
    # imaps1 = (nsplits,Ny,Nx) -- the split maps
    # weights1 = (nsplits,Ny,Nx) -- the weights for the splits

    # select the PS function
    if pfunc is None: pfunc = naive_power
    if not(np.all(np.isfinite(kmaps1))): raise ValueError
    if not(np.all(np.isfinite(weights1))): raise ValueError
    if kmaps2 is not None:
        if not(np.all(np.isfinite(kmaps2))): raise ValueError
        if not(np.all(np.isfinite(weights1))): raise ValueError

    if coadd_estimator:
        # N^{GO} estimator
        assert kmaps1.ndim==3
        nsplits = kmaps1.shape[0]
        if kmaps2 is not None: assert nsplits==kmaps2.shape[0]
        noise = np.sum(pfunc(kmaps1,weights1,kmaps2,weights2),axis=0)
        if not(np.all(np.isfinite(noise))): raise ValueError
        return noise / nsplits / (nsplits-1.)
    else:
        # Cross and auto power
        assert kmaps1.ndim==3
        nsplits = kmaps1.shape[0]
        if kmaps2 is None:
            assert weights2 is None
            kmaps2 = kmaps1
            weights2 = weights1
        else: assert kmaps2.shape[0]==nsplits
        ncrosses = 0.
        nautos = 0.
        crosses = 0.
        autos = 0.
        for i in range(nsplits):
            for j in range(i,nsplits):
                ret = pfunc(kmaps1[i],weights1[i],kmaps2[j],weights2[j])
                if i!=j:
                    crosses = crosses + ret
                    ncrosses += 1
                else:
                    autos = autos + ret
                    nautos += 1
        crosses = crosses / ncrosses
        autos = autos / nautos
        return (autos-crosses)/nsplits


def get_n2d(ffts,wmaps,plot_fname=None,coadd_estimator=False,dtype=None):
    assert np.all(np.isfinite(ffts))
    assert np.all(np.isfinite(wmaps))
    shape,wcs = ffts.shape[-2:],ffts.wcs
    modlmap = enmap.modlmap(shape,wcs)
    Ny,Nx = shape[-2:]
    nfreqs = ffts.shape[0]
    npol = ffts.shape[2]
    ncomps = nfreqs * npol


    def ncomp_to_freq_pol(index):
        ifreq = index // 3
        ipol = index % 3
        return ifreq, ipol

    n2d = enmap.zeros((ncomps,ncomps,Ny,Nx),wcs,dtype=dtype) # WARNING: type)
    pols = ['I','Q','U']
    for i in range(ncomps):
        for j in range(i,ncomps):
            ifreq,ipol = ncomp_to_freq_pol(i)
            jfreq,jpol = ncomp_to_freq_pol(j)

            isplits = ffts[ifreq,:,ipol]
            iwts = wmaps[ifreq,:,0]
            if i!=j:
                jsplits = ffts[jfreq,:,jpol]
                jwts = wmaps[jfreq,:,0]
            else:
                jsplits = None ; jwts = None

            n2d[i,j] = noise_power(isplits,iwts,jsplits,jwts,
                                        coadd_estimator=coadd_estimator,pfunc=naive_power)
            if i!=j: n2d[j,i] = n2d[i,j]
            if plot_fname is not None: 
                plot("%s_%d_%s_%d_%s" \
                     % (plot_fname,ifreq,pols[ipol],
                        jfreq,pols[jpol]),
                     enmap.enmap(np.arcsinh(np.fft.fftshift(n2d[i,j])),wcs),dg=1,quantile=0)

    return n2d

    
    

def smooth_ps(ps,dfact=(16,16),radial_fit_lmin=300,
              radial_fit_lmax=8000,radial_fit_wnoise_annulus=500,
              radial_fit_annulus=20,radial_pairs=[],plot_fname=None,fill_lmax=None,nsplits=None,delta_ell=None,log=False):
    from tilec import covtools
    ncomps = ps.shape[0]
    assert ncomps==ps.shape[1]
    sps = ps*0.
    modlmap = ps.modlmap()
    wcs = ps.wcs

    # Do diagonals first
    for i in range(ncomps):
        if radial_pairs is None: do_radial = True
        else: do_radial = True if (i,i) in radial_pairs else False

        if log:
            dnoise,_,_ = covtools.noise_block_average(ps[i,i].copy(),nsplits,delta_ell,
                                                      lmin=radial_fit_lmin,lmax=radial_fit_lmax,
                                                      wnoise_annulus=radial_fit_wnoise_annulus,
                                                      bin_annulus=radial_fit_annulus,
                                                      radial_fit=do_radial,fill_lmax=fill_lmax,fill_lmax_width=100,log=True)
        else:
            dnoise,_,_ = covtools.noise_average(ps[i,i].copy(),dfact=dfact,
                                                lmin=radial_fit_lmin,lmax=radial_fit_lmax,
                                                wnoise_annulus=radial_fit_wnoise_annulus,
                                                bin_annulus=radial_fit_annulus,
                                                radial_fit=do_radial,fill_lmax=fill_lmax)
        sps[i,i] = dnoise.copy()
        if plot_fname is not None: plot("%s_unsmoothed_%d_%d" \
                                        % (plot_fname,i,i),
                                        enmap.enmap(np.arcsinh(np.fft.fftshift(ps[i,i])),wcs))
        if plot_fname is not None: plot("%s_smoothed_%d_%d" \
                                        % (plot_fname,i,i),enmap.enmap(np.arcsinh(np.fft.fftshift(sps[i,i])),wcs))


    # Do offdiagonals
    for i in range(ncomps):
        for j in range(i+1,ncomps):
            if radial_pairs is None: do_radial = True
            else: do_radial = True if (i,j) in radial_pairs or (j,i) in radial_pairs else False
            afit = ps[i,j].copy()
            if not(do_radial):
                assert i!=j # otherwise this is a trivial operation
                mul = np.sqrt(sps[i,i]*sps[j,j])
            else:
                mul = 1.
            afit = np.nan_to_num(afit / mul)
            afit[...,modlmap<2] = 0.


            if log:
                dnoise,_,_ = covtools.noise_block_average(afit,nsplits,delta_ell,
                                                          lmin=radial_fit_lmin,lmax=radial_fit_lmax,
                                                          wnoise_annulus=radial_fit_wnoise_annulus,
                                                          bin_annulus=radial_fit_annulus,
                                                          radial_fit=do_radial,fill_lmax=fill_lmax,fill_lmax_width=100,log=False)
            else:
                dnoise,_,_ = covtools.noise_average(afit,dfact=dfact,
                                                    lmin=radial_fit_lmin,lmax=radial_fit_lmax,
                                                    wnoise_annulus=radial_fit_wnoise_annulus,
                                                    bin_annulus=radial_fit_annulus,radial_fit=do_radial,fill_lmax=fill_lmax)
            dnoise = dnoise * mul

            sps[i,j] = dnoise.copy()
            if i!=j: sps[j,i] = dnoise.copy()
            if plot_fname is not None: plot("%s_unsmoothed_%d_%d" \
                                            % (plot_fname,i,j),
                                            enmap.enmap(np.arcsinh(np.fft.fftshift(ps[i,j])),wcs))
            if plot_fname is not None: plot("%s_smoothed_%d_%d" \
                                            % (plot_fname,i,j),
                                            enmap.enmap(np.arcsinh(np.fft.fftshift(sps[i,j])),wcs))
            
    sps[...,modlmap<2] = 0.
            
    return sps



def binary_mask(mask,threshold=0.5):
    m = np.abs(mask)
    m[m<threshold] = 0
    m[m>threshold] = 1
    return m


def get_p1ds(p2d,modlmap,bin_edges,dtype=None):
    p1ds = np.zeros((p2d.shape[0],p2d.shape[0],bin_edges.size-1),dtype=dtype) # WARNING: type)
    for i in range(p2d.shape[0]):
        for j in range(p2d.shape[0]):
            p1ds[i,j] = binned_power(p2d[i,j],modlmap,bin_edges)[1]
    cents = (bin_edges[1:] + bin_edges[:-1])/2.
    return cents,p1ds


def compare_ps(cents,p1ds1,p1ds2,plot_fname=None,err=None):
    from orphics import io

    dpi = 300
    ncomps = p1ds1.shape[0]
    if ncomps==1:
        pols = ['I']
    elif ncomps==3:
        pols = ['150-I','150-Q','150-U']
    elif ncomps==6:
        pols = ['150-I','150-Q','150-U','90-I','90-Q','90-U']

    # auto-corrs and cross-freq-II
    k = 0
    pl = io.Plotter(xlabel = "$\\ell$", ylabel = "$D_{\\ell} (\\mu K^2)$",xyscale='linlog',scalefn=lambda x: x**2./2./np.pi)
    for i in range(p1ds1.shape[0]):
        for j in range(p1ds1.shape[0]):
            if not(i==j or (i==0 and j==3)): continue
            polstring = "%s x %s" % (pols[i],pols[j])
            if err is not None:
                   pl.add_err(cents,p1ds1[i,j],
                                yerr=err[i,j],
                                label=polstring,color="C%d" % (k%10),ls="--",alpha=0.5,markersize=2)
            else:
                   pl.add(cents,p1ds1[i,j],
                            label=polstring,color="C%d" % (k%10),ls="--",alpha=0.5,markersize=2)
            pl.add(cents,p1ds2[i,j],color="C%d" % (k%10))
            

            k = k+1
    pl._ax.set_xlim(30,10000)
    pl.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    #pl._ax.set_ylim(5e-6,3e0)
    pl._ax.set_ylim(1e-1,1e5)
    pl.vline(x=500)
    pl.vline(x=300,ls='-.')
    pl.done(plot_fname+"_power.png", dpi=dpi)

    # ratios wrt data of auto-corrs and cross-freq-II
    k = 0
    pl = io.Plotter(xlabel = "$\\ell$", ylabel = "$N_{\\mathrm{sim}} / N_{\\mathrm{data}}$",xyscale='linlin')
    for i in range(p1ds1.shape[0]):
        for j in range(p1ds1.shape[0]):
            if not(i==j or (i==0 and j==3)): continue
            polstring = "%s x %s" % (pols[i],pols[j])
            if err is not None:
                   pl.add_err(cents,p1ds1[i,j]/p1ds2[i,j],yerr=err[i,j]/p1ds2[i,j],label=polstring,color="C%d" % (k%10),alpha=0.7,markersize=3)
            else:
                   pl.add(cents,p1ds1[i,j]/p1ds2[i,j],label=polstring,color="C%d" % (k%10),alpha=0.7,markersize=3)
            k = k+1
    pl._ax.set_xlim(30,10000)
    pl._ax.set_ylim(0.8,1.2)
    pl.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    pl.vline(x=500)
    pl.vline(x=300,ls='-.')
    pl.hline(y=1)
    pl.hline(y=1.05)
    pl.hline(y=0.95)
    pl.done(plot_fname+"_ratio.png", dpi=dpi)

    # cross-pol
    k = 0
    pl = io.Plotter(xlabel = "$\\ell$", ylabel = "$D_{\\ell}  (\\mu K^2)$",xyscale='linlin',scalefn=lambda x: x**2./2./np.pi)
    for i in range(p1ds1.shape[0]):
        for j in range(i+1,p1ds1.shape[0]):
            if ((i==0 and j==3) or (i==3 and j==0)): continue
            polstring = "%s x %s" % (pols[i],pols[j])
            if err is not None:
                   pl.add_err(cents[::2],p1ds1[i,j][::2],
                                yerr=err[i,j][::2],
                                label=polstring,color="C%d" % (k%10),alpha=1,markersize=2)
            else:
                   pl.add(cents[::2],p1ds1[i,j][::2],
                            label=polstring,color="C%d" % (k%10),alpha=1,markersize=2)
            pl.add(cents,p1ds2[i,j],color="C%d" % (k%10),lw=1,alpha=0.5)
            k = k+1
    pl._ax.set_xlim(30,10000)
    pl._ax.set_ylim(-40,40)
    pl.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    pl.hline(y=0)
    pl.vline(x=500)
    pl.vline(x=300,ls='-.')
    pl.done(plot_fname+"_cross_power.png", dpi=dpi)

def plot(fname,imap,dg=4,grid=False,**kwargs):
    from orphics import io

    img = enplot.plot(enmap.downgrade(imap,dg),grid=grid,**kwargs)
    if fname is None: 
        enplot.show(img)
    else: 
        enplot.write(fname,img)
        print(io.bcolors.OKGREEN+"Saved high-res plot to", fname+io.bcolors.ENDC)

def corrcoef(n2d):
    o2d = n2d.copy()*0. + 1.
    for i in range(n2d.shape[0]):
        for j in range(i+1,n2d.shape[0]):
            o2d[i,j] = n2d[i,j] / np.sqrt(n2d[i,i]*n2d[j,j])
            o2d[j,i] = o2d[i,j].copy()
    return o2d



def plot_corrcoeff(cents,c1ds_data,plot_fname):
    from orphics import io
    dpi = 300
    ncomps = c1ds_data.shape[0]
    if ncomps==3:
        pols = ['150-I','150-Q','150-U']
    elif ncomps==6:
        pols = ['150-I','150-Q','150-U','90-I','90-Q','90-U']

    pl = io.Plotter(xlabel = "$\\ell$", ylabel = "$N_{XY}/\\sqrt{N_{XX}N_{YY}}$",xyscale='linlin')
    for i in range(c1ds_data.shape[0]):
        for j in range(i+1,c1ds_data.shape[0]):
            polstring = "%s x %s" % (pols[i],pols[j])
            pl.add(cents,c1ds_data[i,j],label=polstring)
    pl._ax.set_xlim(30,10000)
    pl._ax.set_ylim(-0.3,0.3)
    pl.hline(y=0.05)
    pl.hline(y=-0.05)
    pl.hline(y=0.01,ls='-.')
    pl.hline(y=-0.01,ls='-.')
    pl.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    pl.hline(y=0,ls='-')
    pl.vline(x=500)
    pl.vline(x=300,ls='-.')
    pl.done("%s_compare_corrcoeff.png" % (plot_fname), dpi=dpi)


"""
API

sgen = SimGen()
sgen.get_cmb(season,array,patch)
sgen.get_noise(season,array,patch)
sgen.get_fg(season,array,patch)
sgen.get_sim(season,array,patch)


I set paths:
act_mr3_maps
planck_hybrid_maps
louis_maps

I enforce a base data model spec
dm = DataModel(paths)
dm.get_split()
dm.get_ivar()
dm.get_beam()

# e.g. Planck hybrid
dm.get_split(None,"143","deep56")


"""


