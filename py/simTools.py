import numpy as np, cmblens.flipper.flipperDict as flipperDict, pickle, \
    cmblens.flipper.liteMap as liteMap 


import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
import liteMapPol

from enlib import enmap
import pdb

def getActpolCmbSim(beamfile, coords, iterationNum, cmbDir, cmbSet = 0, \
                    coordsEpsilonArcmin = np.array([[0,0], [0,0]]), \
                    doBeam = True):
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

        

    coordsForEnmap = (coords + coordsEpsilonArcmin / 60.) \
                     * np.pi / 180.

    # pdb.set_trace()

    for x in range(0, 3):
        thisMap = enmap.read_fits(
            cmbDir + \
            "/cmb_set%02d_%05i/fullskyLensedMap_%s_%05d.fits" \
            % (cmbSet, iterationNum, 'TQU'[x], iterationNum), \
            box = coordsForEnmap, \
            wcs_override = wcs )

        # thisMap = enmap.upgrade(thisMap, 2.)
        
        flipperized[x] = thisMap.to_flipper()

    if doBeam:
        output = liteMapPol.simpleBeamConvolution(flipperized[0], \
                                                  flipperized[1], \
                                                  flipperized[2],\
                                                  beamfile )

    return output
    


def getActpolNoiseSim(noiseSeed, patch, noisePsdDir, mask, verbose=True):
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

    # TF.writeFits(noiseDir+'/noise_2dPSD_v2_T_%s.fits'%(patch), overWrite=True)
    # QF.writeFits(noiseDir+'/noise_2dPSD_v2_Q_%s.fits'%(patch), overWrite=True)
    # UF.writeFits(noiseDir+'/noise_2dPSD_v2_U_%s.fits'%(patch), overWrite=True)

    return output


def getActpolSim(iterationNum = 0, region = 'deep5', 
                 season = 's13', \
                 pa = 'pa1', \
                 freqGHz = 150, \
                 patch = None,\
                 coaddDictFile = 'Coadd_s131415.dict', \
                 coaddDictFilePath = '/global/homes/e/engelen/cmblens/inputParams/', \
                 simToolsDictFile = 'simTools.dict',\
                 simToolsDictFilePath = '/global/homes/e/engelen/lensSims/inputParams/',\
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

# def makeActpolNoiseSim(noiseSeed, patch, noisePsdDir, mask)

        return getActpolNoiseSim(noiseSeed = noiseSeed, \
                                  patch = patch, \
                                  noisePsdDir = sDict['noisePsdDir'],
                                  mask = mask,
                                  verbose = verbose)

    elif simType == 'cmb':
# def getActpolCMBSim(beamfile, coords, iterationNum, cmbDir):
        
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
    
# def getActpolCmbSim(iterationNum = 0, region = 'd56', \
#                      season = 's14', \
#                      pa = 'pa1', \
#                      freqGHz = 150, \
#                      patch = None,\
#                      coaddDictFile = 'Coadd_s131415.dict', \
#                      dictFilePath = '/global/homes/e/engelen/cmblens/inputParams/', \
#                      noisePsdDir = '/global/homes/e/engelen/cmblens/maps/dataMaps/actpolDeep/', \
#                      verbose = True)






