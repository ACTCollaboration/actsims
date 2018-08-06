import numpy as np, flipper.flipperDict as flipperDict, pickle
    # flipper.liteMap as liteMap
import astropy.wcs

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)


from enlib import enmap, fft, powspec, resample

import pdb
import os
import scipy.interpolate

def resample_fft_withbeam(d, n, axes=None, doBeam = False, beamData = None, applyWindow = False):
	"""Resample numpy array d via fourier-reshaping. Requires periodic data.
	n indicates the desired output lengths of the axes that are to be
	resampled. By the fault the last len(n) axes are resampled, but this
	can be controlled via the axes argument.

        van Engelen Jan 2018:
        tweaked version which allows to simultaneaously convolve with a beam file and window function.
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
        #AvE hacking
        if applyWindow:
            wy, wx = enmap.calc_window(fd.shape)
            fd *= wy[:,None]**1 * wx[None,:]**1
	# And transform back
	res  = fft.ifft(fd, axes=axes, normalize=True)
	del fd
	res *= norm
	return res if np.issubdtype(d.dtype, np.complexfloating) else res.real


def getActpolCmbSim(beamfileDict,
                    shape, wcs,
                    #box,  #this is what it used to be..
                    iterationNum, cmbDir, freqs, psa,
                    cmbSet = 0, \
                    doBeam = True, pixelFac = 2, applyWindow = True, verbose = True,
                    inscribe = True): #set to True to inscribe inside the shape of the ACTPol map
    #wcs fix per https://phy-wiki.princeton.edu/polwiki/pmwiki.php?n=ACTpolLensingAnalysisTelecon.SomePreliminarySignalSims.  

    shapeFull,wcsFull = enmap.fullsky_geometry(res=1.0*np.pi/180./60.)

    flipperized = [None] * 3

    # coordsForEnmap = np.array([[coords[1,0], coords[0,1]],[ coords[0,0],coords[1,1]]]) \
    #                  * np.pi / 180.

        
    nTQUs = len('TQU')
    firstTime = True

    for tqui in range(0, 3):
        if verbose:
            print 'getActpolCmbSim(): cutting out %s map ' % 'TQU'[tqui]
        lowresMap = enmap.read_fits(cmbDir + \
                                  "fullskyLensedMap_%s_%05d.fits" \
                                  % ('TQU'[tqui], iterationNum), \
                                  box = enmap.box(shape, wcs), \
                                  wcs_override = wcsFull )

        

        for fi, freq in enumerate(freqs):
            beamData = np.loadtxt(os.path.dirname(os.path.abspath(__file__))+"/"+beamfileDict[psa + '_' + freq] ) if doBeam else None

            if verbose:
                print 'getActpolCmbSim(): upsampling by a factor %d for %s.  Beam is %s, pixel window is %s'\
                    % (pixelFac, freq, ('on' if doBeam else 'off'), ('on ' if applyWindow else 'off'))

            upsampled = resample_fft_withbeam(lowresMap, \
                                              shape, \
                                              doBeam = doBeam, 
                                              beamData = beamData, applyWindow = applyWindow)



            if False: #we used to do this check.  But now the output shape should
                #just be equal to the input size
                oshape, owcs = enmap.scale_geometry(lowresMap.shape, lowresMap.wcs, pixelFac )

                if oshape != upsampled.shape:
                    raise ValueError('shape mismatch.  ' + oshape + ' vs. ' \
                                     + 'upsampled.shape')
            if firstTime:
                
                # output = enmap.enmap(np.zeros([len(freqs), 3] \
                #                               +  (list(actpolShape[-2:]) if inscribe else oshape[-2:]) ), owcs)

                output = enmap.enmap(np.zeros([len(freqs), 3] \
                                              +  list( shape[-2:]) )
                                     , wcs)

                firstTime = False


            output[fi, tqui,:, :] = upsampled
        
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

        #first one is for IXI, QxQ, UXU only
        if False:
            print 'loading' + noisePsdDir + '/bigMatrixNoisePsdsCovSqrtDiags_' + psa + '.fits HACKING' 
            covsqrt = enmap.read_fits(noisePsdDir + '/bigMatrixNoisePsdsCovSqrtDiags_' + psa + '.fits' )
        if False:
            print 'loading' + noisePsdDir + '/bigMatrixNoisePsdsCovSqrt_' + psa + '.fits' 
            covsqrt = enmap.read_fits(noisePsdDir + '/bigMatrixNoisePsdsCovSqrt_' + psa + '.fits' )
        if True:
            print 'loading' + noisePsdDir + '/noisePsds_flattened_covSqrt_' + psa + '.fits' 
            covsqrt = enmap.read_fits(noisePsdDir + '/noisePsds_flattened_covSqrt_' + psa + '.fits' )

        
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

        if verbose:
            print 'getActpolNoiseSim(): done '



        return outMaps
    

    else:
        raise ValueError('older ways of getting the noise maps are deprecated')    



def getActpolForegroundSim(beamfileDict ,
                           shape,
                           wcs,
                           iterationNum ,
                           coordsEpsilonArcmin , 
                           doBeam  , applyWindow,
                           foregroundPowerFile,
                           freqs,
                           psa,
                           foregroundSeed, verbose = True):


    from enlib import curvedsky
    # inPowers = np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)),
    #                                    '../data/',
    #                                    foregroundPowerFile))
    

    #This is rescaled ACT foregrounds best-fit for: TT 95x95, 95x150, 150x150, EE 95x95, 95x150, 150x150, TE 95x95, 95x150, 150x150.

    #This means it is expand = 'row' according to the documentation of powspec.compressed_order()
    #We could add polarization to the sims without the argument ncol = 3 -- set to T only for now.

    temperaturePowers \
        = powspec.read_spectrum(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             '../data/',
                                             foregroundPowerFile),
                                ncol = 3,
                                expand = 'row')

    if verbose:
        print 'getActpolForegroundSim(): Getting foreground map '
    outputTT_90_150 = curvedsky.rand_map((2,shape[0], shape[1]),
                                         wcs,
                                         temperaturePowers,
                                         spin = 0,
                                         seed = foregroundSeed)
    outputFreqs = ['f090', 'f150']

    output = enmap.enmap(np.zeros([len(freqs), 3] + list(shape)), wcs)

    for fi, freq in enumerate(freqs):
        if freq in outputFreqs:
            
            output[fi, 'TQU'.index('T'), :, :]  \
                = outputTT_90_150[outputFreqs.index(freq), :, :]



    if doBeam or applyWindow:
        for fi, freq in enumerate(freqs):
            if doBeam:
                if verbose:
                    print 'getActpolForegroundSim(): Convolving foregrounds with beam for frequency ', freq
                    print os.path.dirname(os.path.abspath(__file__))+"/"+beamfileDict[psa + '_' + freq]
                beamData = np.loadtxt(os.path.dirname(os.path.abspath(__file__))+"/"+beamfileDict[psa + '_' + freq] )

                beam2d = scipy.interpolate.interp1d(beamData[:,0], beamData[:,1],
                                                    bounds_error = False, fill_value = 0.)(output.modlmap())
            else:
                beam2d = np.ones(shape)

            if applyWindow:
                print 'getActpolForegroundSim(): Applying pixel window function for frequency ', freq

                wy, wx = enmap.calc_window(shape)

            else :
                wy = np.ones(shape[-2])
                wx = np.ones(shape[-1])

            #beam-convolve all of T, Q, and U at once.
            #also apply the pixel window functions; this line stolen from enmap.apply_window().
            output[fi] = enmap.ifft(enmap.fft(output[fi]) * beam2d * wy[:,None]**1 * wx[None,:]**1).real

    if verbose:
        print 'getActpolForegroundSim(): done '

    return output


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
                 doBeam = True, applyWindow = True):
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


    #unroll psa names (normally stored as a nested  list of lists)
    psaList = [item for sublist in nDict['psaList'] for item in sublist]


    if psa not in psaList:
        raise ValueError('psa %s not found in psaList; options are ' % (psa ), psaList)

    noiseSeed = psaList.index(psa) * 1000000 + iterationNum 

    #load up one sample map, just to get the shape and wcs info.  Do this for "I" at one frequency
    sampleMap = enmap.read_map(os.path.join(os.path.dirname(os.path.abspath(__file__)))+"/"+nDict['dataMapDir'] + 'totalWeightMap' \
                                                        + 'I' + '_' + psa + '_' + psaFreqs[0]  + '_fromenlib.fits') 


    #Note! Foreground seed is the same for every sky patch, season, and frequency!
    #This is because they are used to generate fullsky alm's
    foregroundSeed =  100000000 + iterationNum 

    if simType == 'noise':

        return getActpolNoiseSim(noiseSeed = noiseSeed, \
                                 psa = psa, \
                                 noisePsdDir = os.path.dirname(os.path.abspath(__file__))+"/"+nDict['dataMapDir'],
                                 freqs = psaFreqs, 
                                 verbose = verbose)

    elif simType == 'cmb':
        
        #note dec comes first in the array ordering for enmap, so follow that here
        return getActpolCmbSim(beamfileDict = sDict['beamNames'],
                               # box = enmap.box(sampleMap.shape, sampleMap.wcs),
                               shape = sampleMap.shape, wcs = sampleMap.wcs,
                               iterationNum = iterationNum,
                               cmbDir = os.path.dirname(os.path.abspath(__file__))+"/"+sDict['cmbDir'],
                               freqs = psaFreqs,
                               psa = psa,
                               cmbSet = 0, 
                               doBeam = doBeam, applyWindow = applyWindow,
                               verbose = verbose)

    elif simType == 'foregrounds':

        return getActpolForegroundSim(beamfileDict = sDict['beamNames'],
                                      shape = sampleMap.shape,
                                      wcs = sampleMap.wcs,
                                      iterationNum = iterationNum,
                                      foregroundPowerFile = sDict['foregroundPowerFile'],
                                      coordsEpsilonArcmin = np.array(sDict['coordsEpsilonArcminCMBSim']),\
                                      doBeam = doBeam, applyWindow = applyWindow,
                                      psa = psa,
                                      freqs = psaFreqs,
                                      foregroundSeed = foregroundSeed,
                                      verbose = verbose)




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

