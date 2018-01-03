import numpy, cmblens.flipper.flipperDict as flipperDict, pickle, \
    cmblens.flipper.liteMap as liteMap


def makeActpolNoiseSim(iterationNum = 0, region = 'd56', \
                       season = 's14', \
                       pa = 'pa1', \
                       freqGHz = 150, \
                       patch = None,\
                       inDictFile = 'Coadd_s131415.dict', \
                       dictFilePath = '/global/homes/e/engelen/cmblens/inputParams/', \
                       noisePsdDir = '/global/homes/e/engelen/cmblens/maps/dataMaps/actpolDeep/', \
                       verbose = True):
#FIXME get rid of absolute pathnames
    """ Return a given noise sim  """

    p = flipperDict.flipperDict()
    p.read_from_file(dictFilePath + inDictFile)

    if patch == None:

        freqStr = (freqGHz if type(freqGHz) is str else 'f%03i' % freqGHz)
        patch = '%s-%s-%s-%s' %(season, region, pa, freqStr)


    #unroll patch names (normally stored as a nested  array of arrays)
    patchList = [item for sublist in p['coaddPatchList'] for item in sublist]

    
    if patch not in patchList:
        raise ValueError('patch %s not found in patchList; options are ' % (patch ), patchList)
    mySeed = patchList.index(patch) * 1000000 + iterationNum 

    tNoise = pickle.load(open(noisePsdDir+'noisePowerIAlt_'+patch+'.pkl'))
    qNoise = pickle.load(open(noisePsdDir+'noisePowerQAlt_'+patch+'.pkl'))
    uNoise = pickle.load(open(noisePsdDir+'noisePowerUAlt_'+patch+'.pkl'))

    #FIXME funny naming for now.  Also assumes T weighting for all
    mask = liteMap.liteMapFromFits(noisePsdDir + 'weightMap_T' + patch + '.fits')

    loc = numpy.where(mask.data == 0.)
    loc2 = numpy.where(mask.data < 0.)
    mask.data = numpy.sqrt(mask.data)

    TF = mask.copy()
    QF = mask.copy()
    UF = mask.copy()
    TF.data *= 0.
    QF.data *= 0.
    UF.data *= 0.

    if verbose:
        print 'makeActpolNoiseSim: setting unique RNG seed of %i' % mySeed
    numpy.random.seed(mySeed)
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

    # TF.writeFits(noiseDir+'/noise_2dPSD_v2_T_%s.fits'%(patch), overWrite=True)
    # QF.writeFits(noiseDir+'/noise_2dPSD_v2_Q_%s.fits'%(patch), overWrite=True)
    # UF.writeFits(noiseDir+'/noise_2dPSD_v2_U_%s.fits'%(patch), overWrite=True)





    return [TF, QF, UF]
