from orphics import maps,io,stats
from enlib import bench
from pixell import enmap
import numpy, numpy as np
import os
from tilec import covtools

def load_maps(freqs,dirpath='/home/msyriac/data/act/maps/mr3',spa='s15_deep56_pa3',nsplits=4):
    """
    returns (nfreqs,nsplits,ncomp,Ny,Nx) imap
    returns (nfreqs,nsplits,1,Ny,Nx) wmap
    """
    shape,wcs = enmap.read_map_geometry("%s/%s_%s_nohwp_night_3pass_4way_set%d_map_srcfree.fits" % (dirpath,spa,freqs[0],0))
    Ny,Nx = shape[-2:]
    nfreqs = len(freqs)
    ncomp = 3
    wmaps = np.empty((nfreqs,nsplits,1,Ny,Nx))
    imaps = np.empty((nfreqs,nsplits,ncomp,Ny,Nx))
    coadds = np.empty((nfreqs,ncomp,Ny,Nx))
    cwmaps = np.empty((nfreqs,1,Ny,Nx))
    for ifreq,freq in enumerate(freqs):
        coadds[ifreq] = enmap.read_map("%s/%s_%s_nohwp_night_3pass_4way_coadd_map_srcfree.fits" % (dirpath,spa,freq))
        cwmaps[ifreq] = enmap.read_map("%s/%s_%s_nohwp_night_3pass_4way_coadd_ivar.fits" % (dirpath,spa,freq))[None]
        for isplit in range(nsplits):
            print(spa,freq,isplit)
            imaps[ifreq,isplit] = enmap.read_map("%s/%s_%s_nohwp_night_3pass_4way_set%d_map_srcfree.fits" % (dirpath,spa,freq,isplit))
            wmaps[ifreq,isplit] = enmap.read_map("%s/%s_%s_nohwp_night_3pass_4way_set%d_ivar.fits" % (dirpath,spa,freq,isplit))[None]
    imaps = enmap.enmap(imaps,wcs,copy=False)
    wmaps = enmap.enmap(wmaps,wcs,copy=False)
    coadds = enmap.enmap(coadds,wcs,copy=False)
    cwmaps = enmap.enmap(cwmaps,wcs,copy=False)
    return imaps,wmaps,coadds,cwmaps



spa = 's13_deep6_pa1'
freqs = ["f150"]
region = "deep6"

# spa = 's15_deep56_pa3'
# freqs = ["f090","f150"]
# region = "deep56"

imaps,wmaps,coadds,cwmaps = load_maps(freqs,spa=spa)
shape,wcs = imaps.shape,imaps.wcs
shape = shape[-2:]
fc = maps.FourierCalc(shape,wcs)
Ny,Nx = shape
nfreqs = imaps.shape[0]
assert nfreqs==len(freqs)
ncomp = 3
pscov = maps.SymMat(ncomp*nfreqs,(Ny,Nx))
spscov = maps.SymMat(ncomp*nfreqs,(Ny,Nx))

xlink = enmap.read_map("/home/msyriac/data/act/maps/steve/%s_mask_run_180323_master_apo_w0.fits" % region)
flattened = imaps * np.sqrt(wmaps)
cflattened = coadds * np.sqrt(cwmaps)
del imaps
del wmaps
del coadds, cwmaps
print("FFT")
with bench.show("fft"):
    kmaps = fc.fft(flattened*xlink) / np.sqrt(np.mean(xlink**2.))
    ckmaps = fc.fft(cflattened*xlink) / np.sqrt(np.mean(xlink**2.))
del flattened,xlink,cflattened
fps = []
for ifreq,freq in enumerate(freqs):
    for i in range(3):
        fps.append((ifreq,i))

radial_fit_lmin = 300
radial_fit_lmax = 8000
radial_fit_wnoise_annulus = 500
radial_fit_annulus = 20
dfact = (16,16)
diagnostic_bin_edges = np.arange(300,8000,80)
binner = stats.bin2D(enmap.modlmap(shape,wcs),diagnostic_bin_edges)
pols = ['I','Q','U']
        
for i in range(len(fps)):
    for j in range(i,len(fps)):
        ifreq,ip = fps[i]
        jfreq,jp = fps[j]

        pid = '%s_%s_%s_%s_%s' % (spa,freqs[ifreq],pols[ip],freqs[jfreq],pols[jp])
        do_radial = (ip==jp)
        
        total,crosses,noise = maps.split_calc(kmaps[ifreq,:,ip,...],jsplits=kmaps[jfreq,:,jp,...],icoadd=ckmaps[ifreq,ip,...],jcoadd=ckmaps[jfreq,jp,...],fourier_calc=fc)
        if do_radial:
            dnoise,_,_ = covtools.noise_average(noise,dfact=dfact,
                                                lmin=radial_fit_lmin,
                                                lmax=radial_fit_lmax,
                                                wnoise_annulus=radial_fit_wnoise_annulus,
                                                bin_annulus=radial_fit_annulus,
                                                radial_fit=do_radial)
        else:
            dnoise = enmap.enmap(np.fft.ifftshift(enmap.project(enmap.downgrade(enmap.enmap(np.fft.fftshift(noise),wcs),dfact),shape,wcs)),wcs)
        print(ifreq,jfreq,ip,jp)
        io.plot_img(maps.ftrans(noise),'%s_n2d.png' % pid,aspect='auto')
        io.plot_img(maps.ftrans(dnoise),'%s_down_n2d.png' % pid,aspect='auto')

        pl = io.Plotter(scalefn = lambda x:x**2./np.pi,xlabel='$\\ell$',ylabel='$D_{\\ell}$')
        cents,n1d = binner.bin(noise)
        pl.add(cents,n1d)
        cents,n1d = binner.bin(dnoise)
        pl.add(cents,n1d,ls="--")
        pl.vline(x=300)
        pl.vline(x=500)
        pl.done('%s_n1d.png' % pid)
        

        pscov[i,j] = noise.copy()
        spscov[i,j] = dnoise.copy()

del kmaps,dnoise,noise,crosses,total
enmap.write_map("%s_unsmoothed_pscov.fits" % spa,enmap.enmap(pscov.data,wcs))
del pscov
enmap.write_map("%s_smoothed_pscov.fits" % spa,enmap.enmap(spscov.data,wcs))
