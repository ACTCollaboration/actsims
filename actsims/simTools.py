import numpy as np, flipper.flipperDict as flipperDict, pickle
    # flipper.liteMap as liteMap
import astropy.wcs

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)


from enlib import enmap, fft, powspec, resample, curvedsky

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


# def getActpolForegroundSim(beamFile, coords, iterationNum, cmbDir, cmbSet = 0, beamfile, coords, iterationNum, cmbDir, cmbSet = 0, \
#                            coordsEpsilonArcmin = np.array([[0,0], [0,0]]), \
#                            doBeam = True, pixelFac = 2):

#     shape,wcs = enmap.fullsky_geometry(res=1.0*np.pi/180./60.)
    
#     if coords[1,1] > 180. :
#         coords[1,1] -= 360.
#     if coords[0,1] > 180. :
#         coords[0,1] -= 360.

    

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
        
        beamData = np.loadtxt(beamfile) if doBeam else None

        upsampled = resample_fft_withbeam(thisMap, \
                                 (thisMap.shape[0] * pixelFac, \
                                  thisMap.shape[1] * pixelFac),
                                          doBeam = doBeam, 
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
    


def clPowerFactor(inEnkiMap):
    #get an fft * conj(fft) into C_l units ( assuming no windowing)  by multiplying by this factor.
    return inEnkiMap.area() / (float(inEnkiMap.shape[0] * inEnkiMap.shape[1]))**2


def freqsInPsas(psa, freqsInArraysDict):
    #figure out which frequencies come in a given 'psa' given that this info is not in the dict file -- just listed by 'pa' there.

    arrayList = freqsInArraysDict.keys()
    for array in arrayList:
        if array in psa:
            psaFreqs = freqsInArraysDict[array]
            continue

    return psaFreqs


def getActpolNoiseSim(noiseSeed, psa, noisePsdDir, freqs, verbose = True,
                      useCovSqrt = True,  killFactor = 30., fillValue = 0.):
    #return array of T, Q, U
    #to-do: these are currently using numpy.FFT and are slow; switch to FFTW if installed.

    #Could also have an alternative version of this using enlib tools.

    if useCovSqrt:
        #in this case it was the CovSqrt's that were saved.  This is based on Mat's code in orphics.
        if verbose:
            print 'getActpolNoiseSim(): getting weight maps; assuming I for all'
        
        iqu = 'I' #FOR NOW


        stackOfMaskMaps = [enmap.read_map(noisePsdDir + 'totalWeightMap' \
                                                        + iqu + '_' + psa + '_' + freq  + '_fromenlib.fits') \
                                         for freq in freqs ]
        thisWcs  = stackOfMaskMaps[0].wcs

        maskMaps = enmap.enmap(np.stack(stackOfMaskMaps), thisWcs)

        covsqrt = enmap.read_fits(noisePsdDir + '/bigMatrixNoisePsdsCovSqrt_' + psa + '.fits' )

        if verbose:
            print 'getActpolNoiseSim(): running map_mul to make random phases'

        #get the right normalization
        covsqrt *= np.sqrt(np.prod(covsqrt.shape[-2:]) / enmap.area(covsqrt.shape[-2:], thisWcs ))

        np.random.seed(noiseSeed)
        kmap = enmap.map_mul(covsqrt, enmap.rand_gauss_harm((covsqrt.shape[0], covsqrt.shape[-2:][0], covsqrt.shape[-2:][1]),
                                                            thisWcs))

        #old way:
        # kmapReshape = kmap.reshape((4, kmap.shape[-2:][0], kmap.shape[-2:][1]))
        # outMaps = enmap.ifft(kmapReshape).real
        # kmap /= sqrt(mask)

        if verbose:
            print 'getActpolNoiseSim(): inverse transforming'

        outMaps = enmap.harm2map(kmap, iau = True)
        #now reshape to have shape [nfreqs, 3, Ny, Nx]
        #The "order = 'F' (row vs. column ordering) is due to the ordering that is done
        #in makeNoisePsds.py for the dichroic arrays,
        #namely I90, Q90, U90, I150, Q150, U150.

        outMaps = outMaps.reshape( len(freqs), outMaps.shape[0] / len(freqs),
                                   outMaps.shape[-2], outMaps.shape[-1],
                                   order = 'F')


        for fi, freq in enumerate(freqs):
            #Note each frequency has its own maskmap, so this loop is important
            thisMaskMap = np.squeeze(maskMaps[fi])
            outMaps[fi, :, :, :] /= np.sqrt(thisMaskMap)

            #Loop over T,Q,U.  Couldn't think of clever way to vectorize this part..
            for z in range(outMaps.shape[-3]):
                outMaps[fi, z][thisMaskMap < thisMaskMap[np.where(np.isfinite(thisMaskMap))].max() / killFactor] \
                    = fillValue



        return outMaps
    

    else:
        raise ValueError('older ways of getting the noise maps are deprecated')    



def getActpolForegroundSim(beamfileDict ,
                           shape,
                           wcs,
                           iterationNum ,
                           coordsEpsilonArcmin , 
                           doBeam  ,
                           foregroundPowerFile,
                           freqs,
                           foregroundSeed):


    # inPowers = np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)),
    #                                    '../data/',
    #                                    foregroundPowerFile))
    

    #This is rescaled ACT foregrounds best-fit for: TT 95x95, 95x150, 150x150, EE 95x95, 95x150, 150x150, TE 95x95, 95x150, 150x150.

    temperaturePowers \
        = powspec.read_spectrum(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             '../data/',
                                             foregroundPowerFile),
                                ncol = 3)

    outputTT_90_150 = curvedsky.rand_map((2,shape[0], shape[1]),
                                         wcs,
                                         temperaturePowers,
                                         spin = 0,
                                         seed = foregroundSeed)
    outputFreqs = ['f090', 'f150']

    output = enmap.enmap(np.zeros([len(freqs), 3] + list(sampleMap.shape)))

    # for fi, freq in enumerate(freqs):

    #     output[]


    # output =     
    # if doBeam:
    #     for fi, freq in enumerate(freqs):
    #         beamData = np.loadtxt(beamfileDict[psa + '_' + freq) 

    #     beam2d = scipy.interpolate.interp1d(beamData[:,0], beamData[:,1])(output.modlmap())

    #     output[i] = enmap.ifft(enmap.fft(output[i]) * beam2d)

    return enmap.enmap(np.stack((output, np.zeros(output.shape), np.zeros(output.shape))))


def getActpolSim(iterationNum = 0, patch = 'deep5', 
                 season = 's13', \
                 array = 'pa1', \
                 psa = None,\
                 noiseDictFile = 'templateInputs.dict', \
                 noiseDictFilePath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'../inputParams/'), \
                 signalDictFile = 'signal.dict',\
                 signalDictFilePath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'../inputParams/'),\
                 verbose = True,\
                 simType = 'noise',
                 cmbSet = 0,
                 doBeam = True):
                 #update the last one to True if possible


#FIXME: get rid of absolute pathnames
    """ Return a given noise sim.  
    Provide either 'psa' or all of 'season',  'pa', and 'patch' .

    Will return a stack of enmaps with shape [n_freqs, 3, Ny, Nx],  where the second element has the elements (T, Q, U).  
    n_freqs will be 1 for pa1 and pa2.

"""

    nDict = flipperDict.flipperDict()
    nDict.read_from_file(noiseDictFilePath + '/' + noiseDictFile)

    sDict = flipperDict.flipperDict()
    sDict.read_from_file(signalDictFilePath + '/' + signalDictFile)
    
    if psa == None: #psa stands for patch, season, 
        psa = '%s_%s_%s' %(patch, season, array)


    #Figure out what frequencies correspond to this array, using the function defined above.
    psaFreqs = freqsInPsas(psa, nDict['freqsInArrays'])


    #unroll psa names (normally stored as a nested  array of arrays)
    psaList = [item for sublist in nDict['psaList'] for item in sublist]


    if psa not in psaList:
        raise ValueError('psa %s not found in psaList; options are ' % (psa ), psaList)

    noiseSeed = psaList.index(psa) * 1000000 + iterationNum 

    #load up one sample map, just to get the shape and wcs info.  Do this for "I" at one frequency
    sampleMap = enmap.read_map(nDict['dataMapDir'] + 'totalWeightMap' \
                                                        + 'I' + '_' + psa + '_' + psaFreqs[0]  + '_fromenlib.fits') 


    #Note! Foreground seed is the same for every sky patch, season, and frequency!
    #This is because they are used to generate fullsky alm's
    foregroundSeed =  100000000 + iterationNum 

    if simType == 'noise':

        return getActpolNoiseSim(noiseSeed = noiseSeed, \
                                 psa = psa, \
                                 noisePsdDir = nDict['dataMapDir'],
                                 freqs = psaFreqs, 
                                 verbose = verbose)

    elif simType == 'cmb':
        
        #note dec comes first in the array ordering for enmap, so follow that here
        return getActpolCmbSim(beamfile = sDict['beamNames'][psa],
                               coords = np.array([[mask.y0, mask.x0],[mask.y1, mask.x1]] ),
                               iterationNum = iterationNum,
                               cmbDir = sDict['cmbDir'], cmbSet = 0, \
                               coordsEpsilonArcmin = np.array(sDict['coordsEpsilonArcminCMBSim']),\
                               doBeam = doBeam)

    elif simType == 'foregrounds':

        return getActpolForegroundSim(beamfileDict = sDict['beamNames'],
                               # coords = np.array([[mask.y0, mask.x0],[mask.y1, mask.x1]] ),
                                      shape = sampleMap.shape,
                                      # wcs = mask.wcs.copy(),
                                      wcs = sampleMap.shape,
                                      iterationNum = iterationNum,
                                      foregroundPowerFile = sDict['foregroundPowerFile'],
                                      coordsEpsilonArcmin = np.array(sDict['coordsEpsilonArcminCMBSim']),\
                                      doBeam = doBeam,
                                      freqs = psaFreqs,
                                      foregroundSeed = foregroundSeed)




    else:
        raise ValueError("bad input")
    





#A BUNCH OF LEGACY STUFF -- WHEN WE WERE LOADING NOISE TEMPLATES IN ANOTHER WAY.  MOVED HERE TO KEEP MAIN CODE CLEANER.
    # sqrtmask = mask.copy()
    # sqrtmask.data = np.sqrt(sqrtmask.data)
    # loc = np.where(mask.data == 0.)
    # loc2 = np.where(mask.data < 0.)
        # for iquind , iqu in enumerate(['I', 'Q', 'U']):
        #     noisePsd = enmap.read_map(noisePsdDir + "noisePower" + iqu + 'Alt_'+patch+'_fromenlib.fits')
            
        #     # noise = enmap.rand_gauss(noisePsd.shape,noisePsd.wcs)
        #     realPart = np.random.randn(noisePsd.shape[0], noisePsd.shape[1]) * np.sqrt(noisePsd / clPowerFactor(noisePsd))
        #     imagPart = np.random.randn(noisePsd.shape[0], noisePsd.shape[1]) * np.sqrt(noisePsd / clPowerFactor(noisePsd))

        #     noise = enmap.ifft(realPart + 1j * imagPart, normalize = True).real / (noisePsd.shape[0] * noisePsd.shape[1])**.5
        #     #Simone had this line, turned off for now.
        #     ##NOTE: if I use np.sqrt(hi), with i=0,1,2,3 I can get 4-way noise realization... which is what Steve wants. 
        #     ##noise = enmap.project(noise,ht.shape,ht.wcs)*np.sqrt(ht)
            
        #     # pdb.set_trace()

        #     noise /= sqrtmask.data
        #     noise[loc] = 0
        #     noise[loc2] = 0
        #     output[iquind] = noise.to_flipper()

