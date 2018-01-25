import numpy as np, cmblens.flipper.flipperDict as flipperDict, pickle, \
    cmblens.flipper.liteMap as liteMap


import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)


from enlib import enmap, fft
from enlib import resample
import pdb
import os
import scipy.interpolate

def resample_fft_withbeam(d, n, axes=None, doBeam = False, beamData = None):
	"""Resample numpy array d via fourier-reshaping. Requires periodic data.
	n indicates the desired output lengths of the axes that are to be
	resampled. By the fault the last len(n) axes are resampled, but this
	can be controlled via the axes argument.

        van Engelen Jan 2018:
        tweaked version which allows to simultaneaously convolve with a beam file.
        stolen from enlib.resample.

        """
	d = np.asanyarray(d)
	# Compute output lengths from factors if necessary
	n = np.atleast_1d(n)
	if axes is None: axes = np.arange(-len(n),0)
	else: axes = np.atleast_1d(axes)
	if len(n) == 1: n = np.repeat(n, len(axes))
	else: assert len(n) == len(axes)
	assert len(n) <= d.ndim
	# Nothing to do?
	if np.all(d.shape[-len(n):] == n): return d
	# Use the simple version if we can. It has lower memory overhead
	if d.ndim == 2 and len(n) == 1 and (axes[0] == 1 or axes[0] == -1):
		return resample_fft_simple(d, n[0])
	# Perform the fourier transform

	fd = fft.fft(d, axes=axes)
        if doBeam:
            beam2d = scipy.interpolate.interp1d(beamData[:,0], beamData[:,1])(d.modlmap())
            fd *= beam2d
	# Frequencies are 0 1 2 ... N/2 (-N)/2 (-N)/2+1 .. -1
	# Ex 0* 1 2* -1 for n=4 and 0* 1 2 -2 -1 for n=5
	# To upgrade,   insert (n_new-n_old) zeros after n_old/2
	# To downgrade, remove (n_old-n_new) values after n_new/2
	# The idea is simple, but arbitrary dimensionality makes it
	# complicated.
	norm = 1.0
	for ax, nnew in zip(axes, n):
		ax %= d.ndim
		nold = d.shape[ax]
		dn   = nnew-nold
		if dn > 0:
			padvals = np.zeros(fd.shape[:ax]+(dn,)+fd.shape[ax+1:],fd.dtype)
			spre  = tuple([slice(None)]*ax+[slice(0,nold//2)]+[slice(None)]*(fd.ndim-ax-1))
			spost = tuple([slice(None)]*ax+[slice(nold//2,None)]+[slice(None)]*(fd.ndim-ax-1))
			fd = np.concatenate([fd[spre],padvals,fd[spost]],axis=ax)
		elif dn < 0:
			spre  = tuple([slice(None)]*ax+[slice(0,nnew//2)]+[slice(None)]*(fd.ndim-ax-1))
			spost = tuple([slice(None)]*ax+[slice(nnew//2-dn,None)]+[slice(None)]*(fd.ndim-ax-1))
			fd = np.concatenate([fd[spre],fd[spost]],axis=ax)
		norm *= float(nnew)/nold
	# And transform back
	res  = fft.ifft(fd, axes=axes, normalize=True)
	del fd
	res *= norm
	return res if np.issubdtype(d.dtype, np.complexfloating) else res.real




def getActpolCmbSim(beamfile, coords, iterationNum, cmbDir, cmbSet = 0, \
                    coordsEpsilonArcmin = np.array([[0,0], [0,0]]), \
                    doBeam = True, pixelFac = 2):
    """ coords is in radians as np.array( [ [dec0,ra0], [dec1, ra1 ] ])"""
    #wcs fix per https://phy-wiki.princeton.edu/polwiki/pmwiki.php?n=ACTpolLensingAnalysisTelecon.SomePreliminarySignalSims.  Coords passed in in degrees

    shape,wcs = enmap.fullsky_geometry(res=1.0*np.pi/180./60.)

    flipperized = [None] * 3

    # coordsForEnmap = np.array([[coords[1,0], coords[0,1]],[ coords[0,0],coords[1,1]]]) \
    #                  * np.pi / 180.

    if coords[1,1] > 180. :
        coords[1,1] -= 360.
    if coords[0,1] > 180. :
        coords[0,1] -= 360.

        
    output = [None] * 3
    
    coordsForEnmap = (coords + coordsEpsilonArcmin / 60.) \
                     * np.pi / 180.



    for x in range(0, 3):
        thisMap = enmap.read_fits(
            cmbDir + \
            "/cmb_set%02d_%05i/fullskyLensedMap_%s_%05d.fits" \
            % (cmbSet, iterationNum, 'TQU'[x], iterationNum), \
            box = coordsForEnmap, \
            wcs_override = wcs )

        # thisMap = enmap.upgrade(thisMap, 2.)
        
        beamData = np.loadtxt(beamfile)

        upsampled = resample_fft_withbeam(thisMap, \
                                 (thisMap.shape[0] * pixelFac, \
                                  thisMap.shape[1] * pixelFac),
                                          doBeam = True, 
                                          beamData = beamData)


        oshape, owcs = enmap.scale_geometry(thisMap.shape, thisMap.wcs, pixelFac )

        if oshape != upsampled.shape:
            raise ValueError('shape mismatch.  ' + oshape + ' vs. ' \
                             + 'upsampled.shape')
        
        outputEnmap = enmap.enmap(upsampled, owcs)
        output[x] = outputEnmap.to_flipper()


        # old version (using flipper):
    # if doBeam:
    #     output = liteMapPol.simpleBeamConvolution(flipperized[0], \
    #                                               flipperized[1], \
    #                                               flipperized[2],\
    #                                               beamfile )

    return output
    


def getActpolNoiseSim(noiseSeed, patch, noisePsdDir, mask, verbose=True):
    #return array of T, Q, U
    #to-do: these are currently using numpy.FFT and are slow; switch to FFTW if installed.

    #Could also have an alternative version of this using enlib tools.

    tNoise = pickle.load(open(noisePsdDir+'noisePowerIAlt_'+patch+'.pkl'))
    qNoise = pickle.load(open(noisePsdDir+'noisePowerQAlt_'+patch+'.pkl'))
    uNoise = pickle.load(open(noisePsdDir+'noisePowerUAlt_'+patch+'.pkl'))



    loc = np.where(mask.data == 0.)
    loc2 = np.where(mask.data < 0.)
    mask.data = np.sqrt(mask.data)

    TF = mask.copy()
    QF = mask.copy()
    UF = mask.copy()
    TF.data *= 0.
    QF.data *= 0.
    UF.data *= 0.

    if verbose:
        print 'getActpolNoiseSim: setting unique RNG seed of %i' % noiseSeed
    np.random.seed(noiseSeed)
    TF.fillWithGRFFromTemplate(tNoise,bufferFactor=1)
    QF.fillWithGRFFromTemplate(qNoise,bufferFactor=1)
    UF.fillWithGRFFromTemplate(uNoise,bufferFactor=1)


    TF.data /= mask.data
    QF.data /= mask.data
    UF.data /= mask.data
    TF.data[loc] = 0.
    QF.data[loc] = 0.
    UF.data[loc] = 0.
    TF.data[loc2] = 0.
    QF.data[loc2] = 0.
    UF.data[loc2] = 0.
    output = [TF, QF, UF]

    return output


def getActpolSim(iterationNum = 0, region = 'deep5', 
                 season = 's13', \
                 pa = 'pa1', \
                 freqGHz = 150, \
                 patch = None,\
                 coaddDictFile = 'Coadd_s131415.dict', \
                 coaddDictFilePath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'../inputParams/'), \
                 simToolsDictFile = 'simTools.dict',\
                 simToolsDictFilePath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'../inputParams/'),\
                 verbose = True,\
                 simType = 'noise',
                 cmbSet = 0,
                 doBeam = True):

#FIXME: get rid of absolute pathnames
    """ Return a given noise sim.  
    Provide either 'patch' or all of 'season', 'freqGHz', 'pa', and 'region' .
    For now returns only a liteMap. """

    cDict = flipperDict.flipperDict()
    cDict.read_from_file(coaddDictFilePath + '/' + coaddDictFile)

    sDict = flipperDict.flipperDict()
    sDict.read_from_file(simToolsDictFilePath + '/' + simToolsDictFile)
    


    if patch == None:
        freqStr = (freqGHz if type(freqGHz) is str else 'f%03i' % freqGHz)
        patch = '%s_%s_%s_%s' %(region, season, pa,  freqStr)


    #unroll patch names (normally stored as a nested  array of arrays)
    patchList = [item for sublist in cDict['coaddPatchList'] for item in sublist]


    if patch not in patchList:
        raise ValueError('patch %s not found in patchList; options are ' % (patch ), patchList)
    noiseSeed = patchList.index(patch) * 1000000 + iterationNum 

    foregroundSeed = patchList.index(patch) * 100000000 + iterationNum 

    mask = liteMap.liteMapFromFits(sDict['noisePsdDir'] + 'weightMap_T_' + patch + '.fits')
    if verbose:
        print 'coords are (x0, y0, x1, y1) = (', mask.x0, mask.y0, mask.x1, mask.y1, ')'

    if simType == 'noise':

        return getActpolNoiseSim(noiseSeed = noiseSeed, \
                                  patch = patch, \
                                  noisePsdDir = sDict['noisePsdDir'],
                                  mask = mask,
                                  verbose = verbose)

    elif simType == 'cmb':

        
        #note dec comes first in the array ordering for enmap, so follow that here
        return getActpolCmbSim(beamfile = sDict['beamNames'][patch],
                               coords = np.array([[mask.y0, mask.x0],[mask.y1, mask.x1]] ),
                               iterationNum = iterationNum,
                               cmbDir = sDict['cmbDir'], cmbSet = 0, \
                               coordsEpsilonArcmin = np.array(sDict['coordsEpsilonArcminCMBSim']),\
                               doBeam = doBeam)

    elif simType == 'foregrounds':
        raise ValueError("not yet implemented")

    else:
        raise ValueError("bad input")
    
