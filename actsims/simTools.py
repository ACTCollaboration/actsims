

import numpy as np, flipper.flipperDict as flipperDict, pickle

import astropy.wcs

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)


from pixell import enmap, fft, powspec, curvedsky# , resample

import pdb
import os
import healpy

#indices for RNG seeds
cmbSeedInd = 0
fgSeedInd = 1
phiSeedInd = 2
noiseSeedInd = 3

def getActpolCmbFgSim(beamfileDict,
                      shape, wcs,
                      iterationNum, cmbDir, freqs, psa,
                      cmbSet = 0, \
                      doBeam = True, applyWindow = True, verbose = True, cmbMaptype = 'LensedCMB',
                      foregroundSeed = 0, simType = 'cmb', foregroundPowerFile = None):


    nTQUs = len('TQU')
    firstTime = True

    output   = enmap.empty((len(freqs), nTQUs,)+shape[-2:], wcs)


    if simType == 'cmb':

        filename = cmbDir + "/fullsky%s_alm_set%02d_%05d.fits" % ( cmbMaptype, cmbSet , iterationNum)
        if verbose:
            print 'getActpolCmbFgSim(): loading CMB a_lms from %s' % filename
        import healpy
        almTebFullskyOnecopy = np.complex128(healpy.fitsfunc.read_alm(filename, hdu = (1,2,3)))

        #Now tile the same for all the frequencies requested (i.e. CMB is same at all frequencies).
        #The beam convolution happens below.
        almTebFullsky = np.tile(almTebFullskyOnecopy, (len(freqs), 1, 1))

        if verbose:
            print 'getActpolCmbFgSim(): done'

    elif simType == 'foregrounds':
        outputFreqs = ['f090', 'f150']

        foregroundPowers \
            = powspec.read_spectrum(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                 '../data/',
                                                 foregroundPowerFile),
                                ncol = 3,
                                expand = 'row')


        if verbose:
            print 'getActpolCmbFgSim(): getting foreground a_lms'
        almTFullsky90and150 = curvedsky.rand_alm(foregroundPowers, seed = foregroundSeed)

        almTebFullsky = np.zeros((len(freqs), nTQUs,) + (len(almTFullsky90and150[-1]),), dtype = np.complex128)

        for fi, freq in enumerate(freqs):
            if freq in outputFreqs:
                almTebFullsky[fi, 'TQU'.index('T'), :]  \
                    = almTFullsky90and150[outputFreqs.index(freq), :]

    #Convolve with beam on full sky
    for fi, freq in enumerate(freqs):
        if doBeam:
            beamFile = os.path.dirname(os.path.abspath(__file__))+"/"+beamfileDict[psa + '_' + freq]
            if verbose:
                print 'getActpolCmbFgSim(): applying beam from %s' % beamFile
            beamData = (np.loadtxt(beamFile ))[:,1]
        else:
            if verbose:
                print 'getActpolCmbFgSim(): not convolving with beam'
            beamData = np.repeat(1., almTebFullsky.shape[-1])

        import healpy

        #couldn't quickly figure out how to vectorize this so loop from 0 to 2.
        for tqui in range(nTQUs):
            almTebFullsky[fi, tqui] = healpy.sphtfunc.almxfl(almTebFullsky[fi, tqui].copy(), beamData)

        #These lines stolen from curvedsky.rand_map
        #doing all freqs at once gives error:
        #sharp.pyx in sharp.execute_dp (cython/sharp.c:12118)()
        #ValueError: ndarray is not C-contiguous
        #so loop over all freqs once for now.
        # curvedsky.alm2map(almTebFullsky, output, spin = [0,2],  verbose = True)
        # outputThisfreq   = enmap.empty(( nTQUs,)+shape[-2:], wcs)
        # curvedsky.alm2map(almTebFullsky[fi,:,:], outputThisfreq, spin = [0,2],  verbose = True)
        # output[fi,...] = outputThisfreq


        curvedsky.alm2map(almTebFullsky[fi,:,:], output[fi,:,:,:], spin = [0,2],  verbose = True)



    if applyWindow:
        from enlib import fft

        #The axes along which to FFT
        axes = [-2, -1]
        if verbose:
            print 'getActpolCmbFgSim(): applying pixel window function'

	fd = fft.fft(output, axes = axes)

        wy, wx = enmap.calc_window(fd.shape)

        twoDWindow = wy[:,None]**1 * wx[None,:]**1

        #Careful, this is quietly multiplying an array with shape [N_freq, N_TQU, N_y, N_x] with one of shape [N_y, N_x]
        fd *= twoDWindow
        if verbose:
            print 'getActpolCmbFgSim(): done'
        output = (fft.ifft(fd, axes = axes, normalize = True)).real
        del fd

    return enmap.ndmap(output, wcs)


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
                      useCovSqrt = True,  killFactor = 30., fillValue = 0., noiseDiagsOnly = False):
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

        if noiseDiagsOnly:
            print 'loading' + noisePsdDir + '/noisePsds_flattened_covSqrtDiags_' + psa + '.fits' 
            covsqrt = enmap.read_fits(noisePsdDir + '/noisePsds_flattened_covSqrtDiags_' + psa + '.fits' )
        elif True:
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




def getActpolSim(iterationNum = 0, patch = 'deep5', 
                 season = 's13', \
                 array = 'pa1', \
                 psa = None,\
                 noiseDictFile = 'templateInputsMr3.dict', \
                 noiseDictFilePath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'../inputParams/'), \
                 signalDictFile = 'signal.dict',\
                 signalDictFilePath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'../inputParams/'),\
                 verbose = True,\
                 simType = 'noise',
                 cmbSet = 0,
                 doBeam = True, applyWindow = True, noiseDiagsOnly = False, cmbMaptype = 'LensedCMB'):
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

    # noiseSeed = psaList.index(psa) * 1000000 + iterationNum 
    noiseSeed = (cmbSet, psaList.index(psa), noiseSeedInd, iterationNum)

    #load up one sample map, just to get the shape and wcs info.  Do this for "I" at one frequency
    sampleMap = enmap.read_map(os.path.join(os.path.dirname(os.path.abspath(__file__)))+"/"+nDict['dataMapDir'] + 'totalWeightMap' \
                                                        + 'I' + '_' + psa + '_' + psaFreqs[0]  + '_fromenlib.fits') 


    #Note! Foreground seed is the same for every sky patch, season, and frequency!
    #This is because they are used to generate fullsky alm's
    foregroundSeed = (cmbSet, 0, fgSeedInd, iterationNum)



    if simType == 'noise':

        return getActpolNoiseSim(noiseSeed = noiseSeed, \
                                 psa = psa, \
                                 noisePsdDir = os.path.dirname(os.path.abspath(__file__))+"/"+nDict['dataMapDir'],
                                 freqs = psaFreqs, 
                                 verbose = verbose, noiseDiagsOnly = noiseDiagsOnly)

    elif simType == 'cmb' or simType == 'foregrounds':
        

        return getActpolCmbFgSim(beamfileDict = sDict['beamNames'],
                                 shape = sampleMap.shape, wcs = sampleMap.wcs,
                                 iterationNum = iterationNum,
                                 cmbDir = os.path.dirname(os.path.abspath(__file__))+"/"+sDict['cmbDir'],
                                 freqs = psaFreqs,
                                 psa = psa,
                                 cmbSet = cmbSet, 
                                 doBeam = doBeam, applyWindow = applyWindow,
                                 verbose = verbose, cmbMaptype = cmbMaptype, foregroundSeed = foregroundSeed,
                                 simType = simType,        foregroundPowerFile = sDict['foregroundPowerFile'])


    else:
        raise ValueError("bad input")
    



