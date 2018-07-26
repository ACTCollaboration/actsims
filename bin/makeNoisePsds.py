#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
# import cmblens.flipper.liteMap as liteMap
# import cmblens.flipper.fftTools as fftTools
# import cmblens.flipper.flipperDict as flipperDict
# import flipper.liteMap as liteMap
# import flipper.fftTools as fftTools
import flipper.flipperDict as flipperDict



import numpy, numpy as np
import pickle
import os
# import liteMapPol
# import cmblens.util as util
# import cmblens.config as conf
# import cmblens.mpi as mpi
# import statsTools
from enlib import enmap, array_ops
import mpiTools
import fnmatch
from mpi4py import MPI
import scipy
import sys

import orphics.maps

import gc

def clPowerFactor(inEnkiMap, flipperNorm = False):
    #get an fft * conj(fft) into C_l units assuming no windowing.
    if flipperNorm:
        return inEnkiMap.area() / (float(inEnkiMap.shape[0] * inEnkiMap.shape[1]))**2
    else:
        return 1





def loadSigurdMapInFlipper(filename, i):

    #i = 0 for I, 1 for Q, 2 for U
    if i not in [0,1,2]:
        raise ValueError('bad input')

    #this is how you load just one of [i,q,u]:
    A = enmap.read_map(filename, sel = np.s_[i,  : , : ])

    A_flipperized = A.to_flipper()

    #Why are we doing this step?  because the flipper version that enlib imports might differ from the one that we are loading here
    output = liteMap.liteMapFromDataAndWCS(A_flipperized.data, A_flipperized.wcs)


    return output



def loadSigurdWeightMapInFlipper(filename, i):
    #load ondiagonal piece.
    #i = 0 for I, 1 for Q, 2 for U
    if i not in [0,1,2]:
        raise ValueError('bad input')

    #this is how you load just one of [i,q,u]:
    A = enmap.read_map(filename, sel = np.s_[i, i,  : , : ])

    A_flipperized = A.to_flipper()


    #Why are we doing this step?  because the flipper version that enlib imports might differ from the one that we are loading here
    output = liteMap.liteMapFromDataAndWCS(A_flipperized.data, A_flipperized.wcs)

    return output


    


def appendSubpatchNumbers(rootnames, N):
    output = []
    for rootname in rootnames:

        output.append([rootname + 'patch%03i' % i for i in range(N)])

    return output

        
coaddDict = '../inputParams/' + sys.argv[1]

## Read in Kappa.dict file parameters
p = flipperDict.flipperDict()
p.read_from_file(coaddDict)
dirList = p['dirList']
mapNameList = p['mapNameList']
endNameList = p['endNameList']
endNamePolList = p['endNamePolList']
weightEndNameList = p['weightEndNameList']
noiseEndNameList = p['noiseEndNameList']

dirPolList = p['dirPolList']
mapNamePolList = p['mapNamePolList']
weightEndNamePolList = p['weightEndNamePolList']
simDir     = p['simDir'] 

## Read in Main.dict file parameters 
dataMapDir = p['dataMapDir']
 #allow this to be in its own pkl file



#van Engelen hacking.
# patchList = sum(patchListI, [])
# patchList = patchListI

psaList = [item for sublist in p['psaList'] for item in sublist]




def getEmptyCopy(inObject):
    output = inObject.copy()
    output.data[:] = 0.
    return output

def coaddMapsWithWeights(mapList, weightList, doNull = False):
    # outputMap = mapList[0]
    # outputWeight = mapList[0])
    outputMap = enmap.ndmap(np.zeros(mapList[0].shape), mapList[0].wcs)
    outputWeight = enmap.ndmap(np.zeros(mapList[0].shape), mapList[0].wcs)

    signs = np.ones(len(mapList))
    if doNull:
        if len(mapList) % 2 != 0:
            raise ValueError('bad input, not even length')

        signs[0:len(mapList)/2] *= -1


    for i in range(len(mapList)):
        outputMap += signs[i] * weightList[i] * mapList[i]
        outputWeight += weightList[i]

    outputMap /= outputWeight

    return outputMap, outputWeight





def meanAutoSpec_enlib(mapList,applySlepianTaper=False,nresForSlepian=3.0, window = None, secondMapList = None, secondWindow = None, enlibNorm = True):
    count = 0 
    for i in xrange(len(mapList)):

        # temp_i = enmap.from_flipper(mapList[i])
        # temp_i_prime = enmap.from_flipper(secondMapList[i]) if secondMapList != None else None

        temp_i = mapList[i].copy()
        if secondMapList != None:
            temp_i_prime = secondMapList[i].copy()


        if window is not None :

            # window_enmap = enmap.from_flipper(window)
            print 'meanCrossSpec: doing windowing, doEnlibWay'
            temp_i *= window / np.sqrt(np.mean(window**2))

            if secondMapList is not None:
                

                temp_i_prime *= secondWindow / np.sqrt(np.mean(secondWindow**2))
            else:
                temp_i *= secondWindow / np.sqrt(np.mean(secondWindow**2))

                    
        temp_i = np.nan_to_num(temp_i)
        if temp_i_prime is not None:
            temp_i_prime = np.nan_to_num(temp_i_prime)

        # p2d = (np.conj(enmap.fft(temp_i, normalize = enlibNorm)) \
        #        *  enmap.fft((temp_i_prime if temp_i_prime != None else temp_i), normalize = enlibNorm)).real * clPowerFactor(temp_i, flipperNorm = (not enlibNorm))

        fc = orphics.maps.FourierCalc(temp_i.shape,temp_i.wcs)
        p2d, kmap1, kmap2 =  fc.power2d(np.nan_to_num(temp_i), np.nan_to_num(temp_i_prime))

        if count == 0:
            p2d0 = p2d.copy()
        else:
            p2d0[:] += p2d[:]
        count += 1

    p2d0[:] /= count

    return p2d0



def meanCrossSpec_enlib(mapList,applySlepianTaper=False,nresForSlepian=3.0, window = None, secondMapList = None, secondWindow = None , enlibNorm = True):

    count = 0 

    for i in xrange(len(mapList)):
        for j in xrange(i ):

            temp_i = mapList[i].copy()
            if secondMapList != None:
                temp_j = secondMapList[j].copy()
            else:
                temp_j = mapList[j].copy()
            # temp_i = enmap.from_flipper( mapList[i] )
            # temp_j = enmap.from_flipper( mapList[j] )


            if window is not None :
                # window_enmap = enmap.from_flipper(window)
                print 'meanCrossSpec: doing windowing, doEnlibWay'
                temp_i *= window / np.sqrt(np.mean(window**2))
                temp_i = np.nan_to_num(temp_i)

            if secondWindow is not None:

                temp_j *= secondWindow / np.sqrt(np.mean(secondWindow**2))
                temp_j = np.nan_to_num(temp_j)
                # temp_j[np.where(np.isnan(temp_j[:]))] = 0.


            # p2d = (np.conj(enmap.fft(temp_i, normalize = enlibNorm)) \
            #        *  enmap.fft(temp_j, normalize = enlibNorm)).real * clPowerFactor(temp_i, flipperNorm = (not enlibNorm))

            fc = orphics.maps.FourierCalc(temp_i.shape,temp_i.wcs)
            p2d, kmap1, kmap2 =  fc.power2d(np.nan_to_num(temp_i),np.nan_to_num( temp_j))


            if count == 0:
                p2d0 = p2d.copy()
            else:
                print 'got here'
                p2d0[:] += p2d[:]
            count += 1


    p2d0[:] /= count

    return p2d0#lBin,clBinCrossMean,powerM

                    
def onedl(N):
    return [None] * N




IQUs = ['I', 'Q', 'U'] 
nIQUs = len(IQUs)
nSplits = p['nSplits']
TQUs = ['T', 'Q', 'U'] 

#for reverse compatibility we need to allow for just one value of nSplits
if type(nSplits) is not list:
    nSplits = [nSplits] * len(psaList)


#the weight maps have other tags, for whatever reason.
iqusForWeight = ['', 'q', 'u']
iqusForInputWeights = ['I', 'QQ', 'UU']


justUseIWeights = p['justUseIWeights'] if 'justUseIWeights' in p else False


print 'len(psaList) = ', len(psaList)
piMin, piMax, delta, rank, size = mpiTools.mpiMinMax(MPI.COMM_WORLD, len(psaList))

# piMin, piMax, delta, rank, size = mpiTools.mpiMinMax(MPI.COMM_WORLD, len(psaList))
print 'piMin %i, piMax %i, delta %i, rank %i, size %i' % (piMin, piMax, delta, rank, size)

print 'rank %i will do patches:' %rank, psaList[piMin:piMax]

import itertools

doAll = True

debugFlag = False
if not justUseIWeights:
    raise ValueError('we need to use I weights for now')

#this code does not know the names of the ACTPol arrays, so load up that list here.
arrayList = p['freqsInArrays'].keys()



if doAll:
    for pi, psa in enumerate(psaList):
        firstTime = True

        piPrime = pi #we do not compute correlations between "psa"s -- this primed index is here for future-proofing.
        psaPrime = psa

        for array in arrayList:
            if array in psa:
                psaFreqs = p['freqsInArrays'][array]
                continue


        for index, (iqu, freq)  in enumerate(itertools.product(IQUs, psaFreqs)):

            for indexPrime, (iquPrime, freqPrime) in enumerate(itertools.product(IQUs, psaFreqs)):
                
                if indexPrime > index:
                    continue

                print psa, ': doing ' , freq , 'with ', freqPrime, '; doing ', iqu , 'with ', iquPrime


                mapsForSimgen = onedl(nSplits[pi])
                weightMapsForSimgen = onedl(nSplits[pi])
                actMaps = onedl(nSplits[pi])
                actMapsPrime = onedl(nSplits[pi])
                mapsForPsds = onedl(nSplits[pi])
                mapsForPsdsPrime = onedl(nSplits[pi])
                
                # piPrime = psaList.index()
                mapsForSimgenPrime = onedl(nSplits[pi])

                # if debugFlag and (IQUs.index(iqu) > 0 or IQUs.index(iquPrime) > 0):
                #     continue
                if debugFlag and (IQUs.index(iqu) != IQUs.index(iquPrime)):
                    continue

                for s in range(nSplits[pi]):

                    if p['isEnkiArr'][pi]:

                        mapNameFull = dirList[pi] + mapNameList[pi].format(freq) + str(s) +  endNameList[pi] + '.fits'
                        actMaps[s]  = enmap.read_map(mapNameFull,
                                                          sel = np.s_[IQUs.index(iqu), :, :])
                        # if p['useWeightMap']:
                        weightNameFull = dirList[pi] + mapNameList[pi].format(freq) + str(s) +  weightEndNameList[pi] + '.fits'
                        weightMapsForSimgen[s] = enmap.read_map(weightNameFull)

                        if len(weightMapsForSimgen[s].shape)  > 2:
                        
                            weightMapsForSimgen[s] = enmap.read_map(weightNameFull, sel = np.s_[0, 0, :, :])
                        # else:
                        #     noiseNameFull =  dirList[pi] + mapNameList[pi].format(freq) + str(s) +  noiseEndNameList[pi] + '.fits'
                        #     weightMapsForSimgen[s] = 1./enmap.read_map(noiseNameFull)**2


                        #FLATTEN
                        mapsForSimgen[s] = actMaps[s] * np.sqrt(weightMapsForSimgen[s])
                        mapsForPsds[s] = actMaps[s] * (weightMapsForSimgen[s])
                        
                        ###################################

                        mapNameFullPrime = dirList[piPrime] + mapNameList[piPrime].format(freqPrime) + str(s) +  endNameList[piPrime] + '.fits'
                        actMapsPrime[s] = enmap.read_map(mapNameFullPrime,
                                                          sel = np.s_[IQUs.index(iquPrime), :, :])
                        # if p['useWeightMap']:

                        #note hits maps don't seem to exist for sigurd's maps, so do weight maps no matter what.
                        weightNameFullPrime = dirList[piPrime] + mapNameList[piPrime].format(freqPrime) + str(s)   \
                                                  + weightEndNameList[piPrime] + '.fits'

                        weightMapForSimgenPrime = enmap.read_map(weightNameFullPrime)

                        if len(weightMapForSimgenPrime.shape)  > 2:
                            weightMapForSimgenPrime = enmap.read_map(weightNameFullPrime, sel = np.s_[0, 0, :, :])

                        # else:
                        #     noiseNameFullPrime =  dirList[piPrime] + mapNameList[piPrime].format(freqPrime) + str(s) \
                        #                           + noiseEndNameList[piPrime] + '.fits'
                        #     weightMapsForSimgenPrime[s] = 1./enmap.read_map(noiseNameFullPrime)**2

                        #FLATTEN
                        mapsForSimgenPrime[s] = actMapsPrime[s] * np.sqrt(weightMapForSimgenPrime)
                        mapsForPsds[s] = actMapsPrime[s] * (weightMapForSimgenPrime)

                    # elif p['isPlanckArr'][pi]:

                    #     mapNameFull = dirList[pi] + mapNameList[pi].format(freq) + str(s) +  endNameList[pi] + '.fits'
                    #     actMaps[s]  = enmap.read_map(mapNameFull,
                    #                                       sel = np.s_[IQUs.index(iqu), :, :])
                    #     # if p['useWeightMap']:
                    #     weightNameFull = dirList[pi] + mapNameList[pi].format(freq) + str(s) +  weightEndNameList[pi] + '.fits'
                    #     weightMapsForSimgen[s] = enmap.read_map(weightNameFull, sel = np.s_[0, 0, :, :])
                    #     # else:
                    #     #     noiseNameFull =  dirList[pi] + mapNameList[pi].format(freq) + str(s) +  noiseEndNameList[pi] + '.fits'
                    #     #     weightMapsForSimgen[s] = 1./enmap.read_map(noiseNameFull)**2


                    #     #FLATTEN
                    #     mapsForSimgen[s] = actMaps[s] * np.sqrt(weightMapsForSimgen[s])

                    #     ###################################

                    #     mapNameFullPrime = dirList[piPrime] + mapNameList[piPrime].format(freqPrime) + str(s) +  endNameList[piPrime] + '.fits'
                    #     mapsForSimgenPrime[s] = enmap.read_map(mapNameFullPrime,
                    #                                       sel = np.s_[IQUs.index(iquPrime), :, :])
                    #     # if p['useWeightMap']:

                    #     #note hits maps don't seem to exist for sigurd's maps, so do weight maps no matter what.
                    #     weightNameFullPrime = dirList[piPrime] + mapNameList[piPrime].format(freqPrime) + str(s)   \
                    #                               + weightEndNameList[piPrime] + '.fits'

                    #     weightMapForSimgenPrime = enmap.read_map(weightNameFullPrime, sel = np.s_[0, 0, :, :])

                    #     # else:
                    #     #     noiseNameFullPrime =  dirList[piPrime] + mapNameList[piPrime].format(freqPrime) + str(s) \
                    #     #                           + noiseEndNameList[piPrime] + '.fits'
                    #     #     weightMapsForSimgenPrime[s] = 1./enmap.read_map(noiseNameFullPrime)**2

                    #     #FLATTEN
                    #     mapsForSimgen[s] *= np.sqrt(weightMapForSimgenPrime)


                    else:
                        #ONLY DO THIS FOR SIMONE's MAPS - NO LEADING ZEROS
                        freqFile = 'f90' if freq == 'f090' else freq
                        freqFilePrime = 'f90' if freqPrime == 'f090' else freqPrime


                        mapNameFull = dirList[pi] + mapNameList[pi].format(freqFile) + str(s) + endNameList[pi] + '_' + iqu + '.fits'
                        actMaps[s] = enmap.read_map(mapNameFull)

                        if p['useWeightMap']:
                            weightNameFull = dirList[pi] + mapNameList[pi].format(freqFile) + str(s) +  weightEndNameList[pi] + '.fits'
                            weightMapsForSimgen[s] = enmap.read_map(weightNameFull)
                        else:
                            noiseNameFull = dirList[pi] + mapNameList[pi].format(freqFile) + str(s) +  noiseEndNameList[pi] + '.fits'
                            weightMapsForSimgen[s] = 1./enmap.read_map(noiseNameFull)**2

                        #FLATTEN
                        mapsForSimgen[s] = actMaps[s] * np.sqrt(weightMapsForSimgen[s])

                        mapNameFullPrime = dirList[piPrime] + mapNameList[piPrime].format(freqFilePrime) + str(s) \
                                           + endNameList[piPrime] + '_' + iquPrime + '.fits'

                        mapsForSimgenPrime[s] = enmap.read_map(mapNameFullPrime)
                        if not justUseIWeights:
                            raise ValueError('we need to use I weights for now')

                        if p['useWeightMap']:
                            weightNameFullPrime = dirList[piPrime] + mapNameList[piPrime].format(freqFilePrime) + str(s) \
                                                  +  weightEndNameList[piPrime] + '.fits'

                            weightMapPrime = enmap.read_map(weightNameFullPrime)
                        else:
                            noiseNameFullPrime = dirList[piPrime] + mapNameList[piPrime].format(freqFilePrime) + str(s) \
                                                 +  noiseEndNameList[piPrime] + '.fits'
                            weightMapPrime = 1./enmap.read_map(noiseNameFullPrime)**2
                        mapsForSimgenPrime[s] *= np.sqrt(weightMapPrime)
                        

                if firstTime:
                    psdsForSimgen = enmap.ndmap(np.zeros((len(psaFreqs) * len(IQUs),
                                                           len(psaFreqs)  * len(IQUs),
                                                           mapsForSimgen[s].shape[0],
                                                           mapsForSimgen[s].shape[1])), actMaps[0].wcs)

                    noisePsds = enmap.ndmap(np.zeros(psdsForSimgen.shape), psdsForSimgen.wcs)
                    firstTime = False


                for pt, powerType in enumerate(['flattened', 'unflattened']):
                    crossLinkWindow = enmap.read_fits(p['crossLinkDict'][psa])
                    crossLinkWindowPrime = enmap.read_fits(p['crossLinkDict'][psaPrime])

                    mapListToUse = [mapsForSimgen, actMaps][pt]
                    mapPrimeListToUse = [mapsForSimgenPrime, actMapsPrime][pt]

                    outputPsd = [psdsForSimgen, noisePsds][pt]
                    
                    # #second element is a product of all four splits with the crosslinkwindow
                    # powerWindowToUse = [crossLinkWindow, weightMapsForSimgen * crossLinkWindow ][pt]
                    # powerWindowPrimeToUse = [crossLinkWindowPrime, crossLinkWindowPrime * weightMapPrime][pt]

                    meanCrossPowers_enlib = meanCrossSpec_enlib(mapListToUse, \
                                                                secondMapList = mapPrimeListToUse,
                                                                window = crossLinkWindow,
                                                                secondWindow = crossLinkWindowPrime)

                    meanAutoPowers_enlib = meanAutoSpec_enlib(mapListToUse, \
                                                              secondMapList = mapPrimeListToUse,
                                                              window = crossLinkWindow,
                                                              secondWindow = crossLinkWindowPrime)

                    outputPsd[index, indexPrime, :, :] = meanAutoPowers_enlib - meanCrossPowers_enlib
                    outputPsd[indexPrime, index, :, :] = meanAutoPowers_enlib - meanCrossPowers_enlib

                totalWeightMap = enmap.ndmap(np.sum(weightMapsForSimgen, axis = 0), weightMapsForSimgen[0].wcs)
                totalWeightMap.write(dataMapDir + 'totalWeightMap' + iqu + '_' + psa + '_' + freq  + '_fromenlib.fits')


                nullCoadd, weightMapTotal = coaddMapsWithWeights( actMaps, weightMapsForSimgen, doNull = True)
                enmap.write_fits(dataMapDir + "nullMaps" + iqu + '_' + psa + '_' + freq + '.fits', nullCoadd)

                fullCoadd, weightMapTotal = coaddMapsWithWeights( actMaps, weightMapsForSimgen, doNull = False)
                enmap.write_fits(dataMapDir + "coaddMaps" + iqu + '_' + psa + '_' + freq + '.fits', fullCoadd)

                

        for pt, powerType in enumerate(['flattened', 'unflattened']):
            outputPsd = [psdsForSimgen, noisePsds][pt]

            gc.collect()

            import time
            start = time.time()
            print 'saving ', powerType
            # np.save(dataMapDir + 'bigMatrixNoisePsds_' + psa + '.np', bigMatrixNoisePsds )
            enmap.write_fits(dataMapDir + 'noisePsds_%s_' % powerType + psa + '.fits',  outputPsd)
            print 'done', time.time() - start

            print 'making covsqrt of %s - diagonals' % powerType
            start = time.time()
            covsqrtDiagsOnly = enmap.enmap(np.zeros(outputPsd.shape), outputPsd.wcs)        
            for i in range(len(IQUs)):
                covsqrtDiagsOnly[i::3, i::3, : , :] \
                    = array_ops.eigpow(outputPsd[i::3,i::3, :, :], 0.5, axes = [0,1])
            print 'done', time.time() - start

            print 'writing covsqrt of %s - diagonals' % powerType
            start = time.time()
            enmap.write_fits(dataMapDir + '/noisePsds_%s_covSqrtDiags_' %powerType + psa + '.fits', covsqrtDiagsOnly )
            print 'done', time.time() - start

            print 'making covsqrt of %s' % powerType
            start = time.time()
            covsqrt = array_ops.eigpow(outputPsd, 0.5, axes = [0,1])
            print 'done', time.time() - start

            print 'saving covsqrt of %s' % powerType
            start = time.time()
            enmap.write_fits(dataMapDir + '/noisePsds_%s_covSqrt_' + psa + '.fits', covsqrt )
            print 'done'        , time.time() - start


        
        # bigMatrixNoisePsds

# pickle.dump(bigMatrixNoisePsds, open(resultDir + 'bigMatrixNoisePsds.pkl', 'wb'))

# pickle.dump(enmap.multi_pow(bitMatrixNoisePsds,          0.5),
#             open(resultDir + 'noisePsdsCovSqrt.pkl', 'wb'))






# for s in range(nSplits[pi]):
#             print 'patch %i, loading split %i for map %s' % (pi, s, iqu)
            
#             if p['isEnkiArr'][pi] or p['isEnki']:
#                 actMaps[s] =  loadSigurdMapInFlipper(dir+mapName+ str(s) +  EndName + '.fits', i)
#             else:
#                 actMaps[s] = liteMap.liteMapFromFits(
#                     dir+mapName+ str(s) +EndName+'_'+ (str(i) if (p['isEnki'] ) else iqu)  +'.fits')

#             if p['addInputUsed']:
#                 actMaps_i[s] = liteMap.liteMapFromFits(
#                     dir+mapName+ str(s) +weightEndName+'_input_used_' + iqu + '.fits')
                
#                 actMaps[s].data += actMaps_i[s].data
#             if selectSubmaps:

#                 actMaps[s] =  (actMaps[s]).selectSubMap(mapx0,mapx1,mapy0,mapy1)

#             if not justUseIWeights:
#                 if not (p['isEnki'] or p['isEnkiArr'][pi]):
#                     weightmaps[s]  = liteMap.liteMapFromFits(
#                         dir+mapName+ str(s) +weightEndName+'_weights_'+ iqusForInputWeights[i] + '.fits')
#                     if selectSubmaps:

#                         weightmaps[s] = weightmaps[s].selectSubMap(mapx0,mapx1,mapy0,mapy1)

#                 else:
#                     raise ValueError('non-I weights not coded up for Enki')
#             else:
#                 print 'Using I weights for all'
#                 if p['isEnkiArr'][pi]:
#                     weightmaps[s] = loadSigurdWeightMapInFlipper(dir+mapName + str(s) +  weightEndNameList[pi] + '.fits', 0)

#                 else:
#                     # if p['isEnki'] or p['isEnkiArr'][pi]:
#                     #     weightName = dir + mapName + str(s) + p['weightEndNameList'][pi] + '.fits'
#                     # else:
#                         # weightName = dir+mapName+ str(s) +weightEndName+'_weights_'+i qusForInputWeights[0] + '.fits'
#                         # weightName = dir+mapName+ str(s) +weightEndName+''+ iqusForInputWeights[0] + '.fits'
#                         #hack dec 2017
#                     weightName = dir+mapName+ str(s) +weightEndName+ '.fits'       

#                     weightmaps[s]  = liteMap.liteMapFromFits(weightName)
#                     if selectSubmaps:

#                         weightmaps[s] = weightmaps[s].selectSubMap(mapx0,mapx1,mapy0,mapy1)
                    
#             # if s == 3 and p['doJan2016Split3Tweak'] and patch in ['7ar1' , '7ar2' ]:
#             #     weightmaps[s] = liteMap.liteMapFromFits(dir+mapNamePol+'3'+weightEndNamePol+'_nobad_weights_I.fits')

#             # print 'doing tweak'

#             # if (not justUseIWeights) or (justUseIWeights and i == 0):



#             if p['flipQSign'] and iqu == 'Q':
#                 actMaps[s].data *= -1.

#         if False:
#             actMapCoadd, weightMapTotal = coaddMapsWithWeights( actMaps, weightmaps)


#         # if makeSplits:
#         #     for zeroOne in [0, 1]:
#         #         splitsToWrite[zeroOne], splitsToWriteWeights[zeroOne] = coaddMapsWithWeights(actMaps[2 * zeroOne : 2 * zeroOne + 2 ],
#         #                                                       weightmaps[2 * zeroOne : 2 * zeroOne + 2])

#         #         splitsToWrite[zeroOne].writeFits(
#         #             dataMapDir+'dataCoadd_' + tqu + '_'+patch+'_split%i.fits' % (zeroOne + 1) ,
#         #             overWrite=True)

#         #         splitsToWriteWeights[zeroOne].writeFits(
#         #             dataMapDir+'weightMap_' + tqu + '_'+patch+'_split%i.fits' % (zeroOne + 1) ,
#         #             overWrite=True)

#         # actMapCoadd.writeFits(dataMapDir + 'dataCoadd_' + tqu + '_'+patch+'.fits',overWrite=True)

#         # weightMapTotal.data[ weightMapTotal.data < numpy.max(weightMapTotal.data)/20.] = 0.
#         # weightMapTotal.writeFits(dataMapDir + 'weightMap_' +  tqu + '_' + patch + '.fits',\
#         #                              overWrite = True)

#         # if i == 0:
#         #     ell = numpy.arange(20000.)

#         #     #then create a taper window for the power spectrum estimation stuff below.
#         #     weightMapTotalTemperature = weightMapTotal.copy()
#         #     weightMapTotalTemperature = weightMapTotalTemperature.filterFromList([ell,numpy.exp(-(ell/1200.)**2./2.)])
#         #     taper = liteMapPol.initializeCosineWindow(weightMapTotalTemperature,p['taperWidth'],0)
#         #     weightMapTotalTemperature.data *= taper.data
#         #     weightMapTotalTemperature.writeFits(dataMapDir+'mask_'+patch+'.fits',overWrite=True)
#         #     sqrtWeightMapTotal = weightMapTotalTemperature.copy()
#         #     sqrtWeightMapTotal.data = numpy.nan_to_num(numpy.sqrt(sqrtWeightMapTotal.data))
            

#         # if makeSplits:
#         #     noiseMap = splitsToWrite[0].copy()
#         #     noiseMap.data -= splitsToWrite[1].data
#         #     noiseMap.data /= 2.
#         #     noiseMap.writeFits(dataMapDir + 'noiseMap_' + tqu + '_'+patch+'.fits',overWrite=True)

#         if True:

#             if p['useCrossLinkMaps']:
#                 #the old way of doing this
#                 if False:
#                     possibleCrossLinkFilenames = fnmatch.filter(os.listdir(p['crossLinkRoot']),\
#                                                                 patch + '*fits')


#                     if len(possibleCrossLinkFilenames) != 1:
#                         raise ValueError('more or less than 1 possible filename match for crosslink')

#                     powerWindow = liteMap.liteMapFromFits(p['crossLinkRoot'] + possibleCrossLinkFilenames[0])
#                     if selectSubmaps:

#                         powerWindow =  powerWindow.selectSubMap(mapx0,mapx1,mapy0,mapy1)
                
#                 if True:
                    
#                     crossLinkMap = p['crossLinkDict'][patch]
#                     if crossLinkMap != None:  #None is the case where they are not made yet
#                         powerWindow = liteMap.liteMapFromFits(p['crossLinkDict'][patch])
#                         if selectSubmaps:
#                             powerWindow =  powerWindow.selectSubMap(mapx0,mapx1,mapy0,mapy1)
                        

#                     else:
#                         powerWindow = actMaps[s].copy()
#                         powerWindow.data[:] = 1.

#                 # powerWindow.writeFits(dataMapDir + 'crosslinkMap_T_' + patch+'.fits',\
#                 #                       overWrite=True)
                

#             else:
#                 powerWindow = None



#             # mapsForSimgen = statsTools.onedl(nSplits[pi])
#             # for s in range(nSplits[pi]):
#             #     mapsForSimgen[s] = actMaps[s].copy()
#             #     mapsForSimgen[s].data *= sqrtWeightMapTotal.data

#             doFlipperCrosses = False
#             if doFlipperCrosses:
#                 meanAutoPowers = meanAutoSpec(mapsForSimgen, \
#                                               window = powerWindow)
#                 meanCrossPowers = meanCrossSpec(mapsForSimgen, \
#                                                 window = powerWindow)

#                 meanAutoPowers.powerMap -= meanCrossPowers.powerMap
#                 w2 = 0.25#1.#/numpy.mean(wCAlt.data**2.)
#                 meanAutoPowers.powerMap *= w2

#                 if not np.all(meanAutoPowers.powerMap > 0):
#                     raise ValueError('power spectrum is negative somewhere')

#                 print 'writing to ' , dataMapDir+'noisePower' + iqu + 'Alt_'+patch+'.pkl'
#                 pickle.dump(meanAutoPowers, open(dataMapDir+'noisePower' + iqu + 'Alt_'+patch+'.pkl','w'))
#                 print 'done'


#             # meanAutoPowers = meanAutoSpec(mapsForSimgen, \
#             #                               window = powerWindow, doEnlibWay = True)


#             # if estimateCrossesBetweenPatches:
#             #     meanAutoPowersExtended[
#             #     ][] = actMaps[s].copy()







#             mapsForPowerEst = statsTools.onedl(nSplits[pi])
#             for s in range(nSplits[pi]):
#                 mapsForPowerEst[s] = actMaps[s].copy()
#                 mapsForPowerEst[s].data *= weightMapTotalTemperature.data

#             meanAutoPowers = meanAutoSpec(mapsForSimgen)
#             meanCrossPowers = meanCrossSpec(mapsForSimgen)

#             meanAutoPowers.powerMap -= meanCrossPowers.powerMap

#             w2 = 0.25 / numpy.mean(weightMapTotalTemperature.data**2.)
#             meanAutoPowers.powerMap *= w2

#             print 'writing to ', dataMapDir+'dataNoise_' + patch + '_' + iqu + '.pkl'
#             pickle.dump(meanAutoPowers, open(dataMapDir+'dataNoise_' + patch + '_' + iqu + '.pkl','w'))
#             print 'done'







#         ##garbage code left here for posterity
#                 #         if p['doSmoothingOfSimPower']:
#                 # print 'smoothing power with width of ', p['smoothingWidthInPixels'], 'pixels'
#                 # meanAutoPowers.powerMap = scipy.ndimage.filters.gaussian_filter(meanAutoPowers.powerMap, \
#                 #                                                                 p['smoothingWidthInPixels'])




#USED TO SAVE EACH MAP SEPARATELY


                # (meanAutoPowers_enlib - meanCrossPowers_enlib).write(dataMapDir + "noisePower" + iqu + '_' + iquPrime + 'Alt_'\
                #                                                      +psa + '_' + freq + '__' + psaPrime + '_' + freqPrime + '_fromenlib.fits')
                # (meanAutoPowers_enlib).write(dataMapDir + "splitAutoPower" + iqu + '_' + iquPrime + 'Alt_'\
                #                                                      +psa + '_' + freq + '__' + psaPrime + '_' + freqPrime + '_fromenlib.fits')
                # (meanCrossPowers_enlib).write(dataMapDir + "splitCrossPower" + iqu + '_' + iquPrime + 'Alt_'\
                #                                                      +psa + '_' + freq + '__' + psaPrime + '_' + freqPrime + '_fromenlib.fits')

