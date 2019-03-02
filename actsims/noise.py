import numpy as np
import os,sys
from pixell import enmap,enplot,fft as pfft
from soapack import interfaces as sints
from actsims import utils
from orphics import io
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
        self._dm = sints.models[model](region=extract_region)
        self.verbose = verbose

    def save_covsqrt(self,covsqrt,season=None,patch=None,array=None,coadd=True,mask_patch=None):
        pout,cout,sout = get_save_paths(self._model,self._version,
                                        coadd=coadd,season=season,patch=patch,array=array,
                                        overwrite=False,mask_patch=mask_patch)
        fpath = "%s_covsqrt.fits" % (cout)
        enmap.write_map(fpath ,covsqrt)
        print(fpath,covsqrt.shape)

    def load_covsqrt(self,season=None,patch=None,array=None,coadd=True,mask_patch=None):
        pout,cout,sout = get_save_paths(self._model,self._version,coadd=coadd,
                                        season=season,patch=patch,array=array,
                                        overwrite=False,mask_patch=mask_patch)
        fpath = "%s_covsqrt.fits" % (cout)
        ikey = '_'.join([str(x) for x in [season,patch,self._dm.array_freqs[array]]])
        try:
            covsqrt = self._ccache[fpath]
            ivars = self._icache[ikey]
            if self.verbose: print("Loaded cached covsqrt and ivars.")
        except:
            if self.verbose: print("Couldn't find covsqrt and ivars in cache. Reading from disk...")
            ivars = self._dm.get_splits_ivar(season=season,patch=patch,arrays=self._dm.array_freqs[array])
            covsqrt = enmap.read_map(fpath)
            print(fpath,covsqrt.shape)
            if len(self._ccache.keys())<self.ncache: 
                self._ccache[fpath] = covsqrt
                self._icache[ikey] = ivars
        return covsqrt,ivars
        
    def generate_sim(self,season=None,patch=None,array=None,seed=None,mask_patch=None,binary_percentile=10.):
        covsqrt,ivars = self.load_covsqrt(season=season,patch=patch,array=array,mask_patch=mask_patch)
        sims = generate_noise_sim(covsqrt,ivars,seed=seed,binary_percentile=binary_percentile)
        return sims


    def save_sims(self,sims,season,patch,array,mask_patch,coadd=True):
        pout,cout,sout = get_save_paths(self._model,self._version,coadd=coadd,
                                        season=season,patch=patch,array=array,
                                        overwrite=False,mask_patch=mask_patch)
        insplits = self._dm.get_nsplits(season,patch,array)
        freqs = self._dm.array_freqs[array]
        nfreqs,nsplits,npol,Ny,Nx = sims.shape
        assert nsplits == insplits
        assert len(freqs) == nfreqs
        for i in range(nfreqs):
            iarray = freqs[i]
            for j in range(nsplits):
                fname = sout+os.path.basename(self._dm.get_split_fname(season,patch,iarray,j,srcfree=True))
                enmap.write_map(fname,sims[i,j,:,:,:])
        


def get_save_paths(model,version,coadd,season=None,patch=None,array=None,mkdir=False,overwrite=False,mask_patch=None):
    paths = sints.dconfig['actsims']

    assert paths['plot_path'] is not None
    assert paths['covsqrt_path'] is not None
    assert paths['trial_sim_path'] is not None

    # Prepare output dirs
    pdir = "%s/%s/" % (paths['plot_path'] ,version) 
    cdir = "%s/%s/" % (paths['covsqrt_path'] ,version)
    sdir = "%s/%s/" % (paths['trial_sim_path'] ,version)
    
    if mkdir:
        exists1 = utils.mkdir(pdir)
        exists2 = utils.mkdir(cdir)
        exists3 = utils.mkdir(sdir)
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


def get_n2d_data(splits,ivars,mask_a,coadd_estimator=False,flattened=False,plot_fname=None):
    assert np.all(np.isfinite(splits))
    assert np.all(np.isfinite(ivars))
    assert np.all(np.isfinite(mask_a))

    if coadd_estimator:
        coadd,_ = get_coadd(splits,ivars,axis=1)
        data  = splits - coadd[:,None,...]
        del coadd
    else:
        data = splits

    assert np.all(np.isfinite(data))
    if flattened:
        ffts = enmap.fft(data*mask_a*np.sqrt(ivars),normalize="phys")
        if plot_fname is not None: plot(plot_fname+"_fft_maps",data*mask_a*ivars)
        wmaps = mask_a + enmap.zeros(ffts.shape,mask_a.wcs,dtype=sints.dtype) # WARNING: type
        del ivars, data, splits
    else:
        assert np.all(np.isfinite(data*mask_a*ivars))
        ffts = enmap.fft(data*mask_a*ivars,normalize="phys")
        if plot_fname is not None: plot(plot_fname+"_fft_maps",data*mask_a*ivars)
        wmaps = ivars * mask_a
        del ivars, data, splits
    n2d = get_n2d(ffts,wmaps,coadd_estimator=coadd_estimator,plot_fname=plot_fname)
    assert np.all(np.isfinite(n2d))
    return n2d


def generate_noise_sim(covsqrt,ivars,binary_percentile=10.,seed=None):
    if isinstance(seed,int): seed = [seed]
    assert np.all(np.isfinite(covsqrt))
    eshape,ewcs = ivars.shape,ivars.wcs
    shape,wcs = covsqrt.shape,covsqrt.wcs
    Ny,Nx = shape[-2:]
    ncomps = covsqrt.shape[0]
    assert ncomps==covsqrt.shape[1]
    assert ncomps % 3 == 0
    nfreqs = ncomps // 3
    wmaps = ivars
    nsplits = wmaps.shape[1]

    if sints.dtype is np.float32: ctype = np.complex64 
    elif sints.dtype is np.float64: ctype = np.complex128 

    # Old way with loop
    kmap = []
    for i in range(nsplits):
        if seed is None:
            np.random.seed(None)
        else:
            np.random.seed(seed+[i])
        rmap = enmap.rand_gauss_harm((ncomps, Ny, Nx),covsqrt.wcs).astype(ctype)
        kmap.append( enmap.map_mul(covsqrt, rmap) )
    del covsqrt, rmap
    kmap = enmap.enmap(np.stack(kmap),wcs)
    outmaps = enmap.extract(enmap.ifft(kmap, normalize="phys").real,eshape,wcs)
    del kmap

    # Need to test this more ; it's only marginally faster and has different seed behaviour
    # covsqrt = icovsqrt 
    # np.random.seed(seed)
    # rmap = enmap.rand_gauss_harm((nsplits,ncomps,Ny, Nx),covsqrt.wcs)
    # kmap = enmap.samewcs(np.einsum("abyx,cbyx->cayx", covsqrt, rmap),rmap)
    # outmaps = enmap.ifft(kmap, normalize="phys").real

    assert np.all(np.isfinite(outmaps))
    # Divide by hits
    for ifreq in range(nfreqs):
        outmaps[:,ifreq*3:(ifreq+1)*3,...] = outmaps[:,ifreq*3:(ifreq+1)*3,...] / np.sqrt(wmaps[ifreq,...]) *np.sqrt(nsplits)

    if binary_percentile is not None:
        # Sanitize by thresholding and binary masking
        for ifreq in range(nfreqs):
            for isplit in range(nsplits):
                win = wmaps[ifreq,isplit,0,...]
                #bmask = binary_mask(win,threshold = np.percentile(win,binary_percentile))
                #outmaps[isplit,ifreq*3:(ifreq+1)*3,bmask==0] = 0
                outmaps[isplit,ifreq*3:(ifreq+1)*3,win==0.] = 0 # FIXME: make 0 comparison robust

    retmaps = outmaps.reshape((nsplits,nfreqs,3,Ny,Nx)).swapaxes(0,1)
    assert np.all(np.isfinite(retmaps))
    return retmaps

    
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


def get_n2d(ffts,wmaps,plot_fname=None,coadd_estimator=False):
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

    n2d = enmap.zeros((ncomps,ncomps,Ny,Nx),wcs,dtype=sints.dtype) # WARNING: type)
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
                     enmap.enmap(np.arcsinh(np.fft.fftshift(n2d[i,j])),wcs))

    return n2d

    
    

def smooth_ps(ps,dfact=(16,16),radial_fit_lmin=300,
              radial_fit_lmax=8000,radial_fit_wnoise_annulus=500,
              radial_fit_annulus=20,radial_pairs=[],plot_fname=None):
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
        dnoise,_,_ = covtools.noise_average(ps[i,i].copy(),dfact=dfact,
                                            lmin=radial_fit_lmin,lmax=radial_fit_lmax,
                                            wnoise_annulus=radial_fit_wnoise_annulus,
                                            bin_annulus=radial_fit_annulus,radial_fit=do_radial)
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
            dnoise,_,_ = covtools.noise_average(afit,dfact=dfact,
                                                lmin=radial_fit_lmin,lmax=radial_fit_lmax,
                                                wnoise_annulus=radial_fit_wnoise_annulus,
                                                bin_annulus=radial_fit_annulus,radial_fit=do_radial)
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


def get_p1ds(p2d,modlmap,bin_edges):
    p1ds = np.zeros((p2d.shape[0],p2d.shape[0],bin_edges.size-1),dtype=sints.dtype) # WARNING: type)
    for i in range(p2d.shape[0]):
        for j in range(p2d.shape[0]):
            p1ds[i,j] = binned_power(p2d[i,j],modlmap,bin_edges)[1]
    cents = (bin_edges[1:] + bin_edges[:-1])/2.
    return cents,p1ds


def compare_ps(cents,p1ds1,p1ds2,plot_fname=None,err=None):

    dpi = 300
    ncomps = p1ds1.shape[0]
    if ncomps==3:
        pols = ['150-I','150-Q','150-U']
    elif ncomps==6:
        pols = ['90-I','90-Q','90-U','150-I','150-Q','150-U']

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
    pl.done(plot_fname+"_cross_power.png", dpi=dpi)

def plot(fname,imap,dg=4,grid=False,**kwargs):
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
    dpi = 300
    ncomps = c1ds_data.shape[0]
    if ncomps==3:
        pols = ['150-I','150-Q','150-U']
    elif ncomps==6:
        pols = ['90-I','90-Q','90-U','150-I','150-Q','150-U']

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


