#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
# import cmblens.flipper.liteMap as liteMap
# import cmblens.flipper.fftTools as fftTools
# import cmblens.flipper.flipperDict as flipperDict
# import flipper.liteMap as liteMap
# import flipper.fftTools as fftTools
from actsims import flipperDict



import numpy, numpy as np
import pickle
import os
# import liteMapPol
# import cmblens.util as util
# import cmblens.config as conf
# import cmblens.mpi as mpi
# import statsTools
from enlib import enmap, array_ops
from actsims import mpiTools
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

        signs[0:len(mapList)//2] *= -1


    for i in range(len(mapList)):
        outputMap += signs[i] * weightList[i] * mapList[i]
        outputWeight += weightList[i]

    outputMap /= outputWeight

    return outputMap, outputWeight





def meanAutoSpec_enlib(mapList,applySlepianTaper=False,nresForSlepian=3.0, window = None, secondMapList = None, secondWindow = None, enlibNorm = True):
    count = 0 
    for i in range(len(mapList)):

        # temp_i = enmap.from_flipper(mapList[i])
        # temp_i_prime = enmap.from_flipper(secondMapList[i]) if secondMapList != None else None

        temp_i = mapList[i].copy()
        if secondMapList != None:
            temp_i_prime = secondMapList[i].copy()


        if window is not None :

            # window_enmap = enmap.from_flipper(window)
            print('meanCrossSpec: doing windowing, doEnlibWay')
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

    for i in range(len(mapList)):
        for j in range(i ):

            temp_i = mapList[i].copy()
            if secondMapList != None:
                temp_j = secondMapList[j].copy()
            else:
                temp_j = mapList[j].copy()
            # temp_i = enmap.from_flipper( mapList[i] )
            # temp_j = enmap.from_flipper( mapList[j] )


            if window is not None :
                # window_enmap = enmap.from_flipper(window)
                print('meanCrossSpec: doing windowing, doEnlibWay')
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
                print('got here')
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


print('len(psaList) = ', len(psaList))
piMin, piMax, delta, rank, size = mpiTools.mpiMinMax(MPI.COMM_WORLD, len(psaList))

# piMin, piMax, delta, rank, size = mpiTools.mpiMinMax(MPI.COMM_WORLD, len(psaList))
print('piMin %i, piMax %i, delta %i, rank %i, size %i' % (piMin, piMax, delta, rank, size))
print('rank %i will do patches:' %rank, psaList[piMin:piMax])

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

                print(pi,psa,index,indexPrime,iqu,iquPrime,freq,freqPrime)
                continue

                print(psa, ': doing ' , freq , 'with ', freqPrime, '; doing ', iqu , 'with ', iquPrime)


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
                        mapsForPsds[s] = actMaps[s] * weightMapsForSimgen[s] / np.mean(weightMapsForSimgen[s]**2)
                        
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
                        mapsForPsdsPrime[s] = actMapsPrime[s] * (weightMapForSimgenPrime) / np.mean(weightMapForSimgenPrime**2)
                        

                if firstTime:
                    psdsForSimgen = enmap.ndmap(np.zeros((len(psaFreqs) * len(IQUs),
                                                           len(psaFreqs)  * len(IQUs),
                                                           mapsForSimgen[s].shape[0],
                                                           mapsForSimgen[s].shape[1])), actMaps[0].wcs)

                    noisePsds = enmap.ndmap(np.zeros(psdsForSimgen.shape), psdsForSimgen.wcs)
                    firstTime = False


                for pt, powerType in enumerate(['flattened']):
                    #, 'unflattened']):
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

                

        for pt, powerType in enumerate(['flattened']):
        #, 'unflattened']):
            outputPsd = [psdsForSimgen, noisePsds][pt]

            gc.collect()

            import time
            start = time.time()
            print('saving ', powerType)
            # np.save(dataMapDir + 'bigMatrixNoisePsds_' + psa + '.np', bigMatrixNoisePsds )
            enmap.write_fits(dataMapDir + 'noisePsds_%s_' % powerType + psa + '.fits',  outputPsd)
            print('done', time.time() - start)


            print('making covsqrt of %s - diagonals' % powerType)
            start = time.time()
            covsqrtDiagsOnly = enmap.enmap(np.zeros(outputPsd.shape), outputPsd.wcs)        
            for i in range(len(IQUs)):
                covsqrtDiagsOnly[i::3, i::3, : , :] \
                    = array_ops.eigpow(outputPsd[i::3,i::3, :, :], 0.5, axes = [0,1])
            print('done', time.time() - start)

            print('writing covsqrt of %s - diagonals' % powerType)
            start = time.time()
            enmap.write_fits(dataMapDir + '/noisePsds_%s_covSqrtDiags_' %powerType + psa + '.fits', covsqrtDiagsOnly )
            print('done', time.time() - start)

            print('making covsqrt of %s' % powerType)
            start = time.time()
            covsqrt = array_ops.eigpow(outputPsd, 0.5, axes = [0,1])
            print('done', time.time() - start)

            print('saving covsqrt of %s' % powerType)
            start = time.time()
            enmap.write_fits(dataMapDir + '/noisePsds_%s_covSqrt_' %powerType + psa + '.fits', covsqrt )
            print('done'        , time.time() - start)


        
