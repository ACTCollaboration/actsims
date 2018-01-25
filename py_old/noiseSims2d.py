# #!/usr/bin/env python  
# import sys
# sys.path.append('../tools')
# sys.path.append('../inputParams')

# from flipper import *
# import pyfits
# import pickle
# import matplotlib.pyplot as plt
# import healpy as H
# import math
# import liteMapPol
# import time
# import mapTools
# import sys
# import inpaintTools


# ## Define parameters	  
# TCMB = 2.726e6

# ## If no .dict file are given as input, read in PrepareMaps.dict as default
# if len(sys.argv) < 2:
#     mapDict = '../inputParams/PrepareMaps.dict'
# else:
#     mapDict = '../inputParams/' + sys.argv[1]  # put desired .dict file in inputParams dir    
# if len(sys.argv) < 3:
#     mainDict = '../inputParams/Main.dict'
# else:
#     mainDict = '../inputParams/' + sys.argv[2]

# ## Read in PrepareMaps.dict file parameters 
# p = flipperDict.flipperDict()
# p.read_from_file(mapDict)
# window = p['window']
# applyWeights = p['applyWeights']
# applyCalib = p['applyCalib']
# removePtSources = p['removePtSources']
# addNoise = p['addNoise']
# addForegrounds = p['addForegrounds']
# inPaintClusters = p['inPaintClusters']
# eps = p['inPaintEpsilon']


# ## Read in Main.dict file parameters															   
# p.read_from_file(mainDict)
# data = p['data']
# resultDir = p['resultDir']
# dataMapDir = p['dataMapDir']
# patchList = p['patchList']

# # patchList = ['0']
# start = time.time()



# from mpi4py import MPI

# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# size = comm.Get_size()

# print rank, size


# iStart = 0
# iStop = 1

# delta = (iStop - iStart) / size
# if delta == 0:
# 	raise ValueError, 'Too many processors for too small a  loop!'

# print delta
# iMin = iStart+rank*delta
# iMax = iStart+(rank+1)*delta
# if iMax>iStop:
# 	iMax = iStop
# elif (iMax > (iStop - delta)) and iMax <iStop:
# 	iMax = iStop
from flipper import *
import pickle


print "Reading dict file"
p = flipperDict.flipperDict()
p.read_from_file(sys.argv[1])

from mpi4py import MPI

#MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
    
print rank, size

# print 'temporary 10'
# iMax = 10

iStart = p['iStart']
iStop = p['iStop']


delta = (iStop - iStart)/size

if delta == 0:
	raise ValueError, 'Too many processors for too small a  loop!'

print delta
iMin = iStart+rank*delta
iMax = iStart+(rank+1)*delta



numCMBsets = p['numCMBsets']

noiseDirRoot = p['noiseDirRoot']


# noisePsdDir = '/home/r/rbond/bsherwin/dis2/septemberLensingVersion/lensRecon/maps/dataMaps/actpolDeep/'
# noisePsdDir = '/scratch2/r/rbond/engelen/cmblens/maps/dataMaps/actpolDeep/'
noisePsdDir = p['noisePsdDir']



weightMapDir = noisePsdDir

# patchList = ['5s1ar1', '6s1ar1', '6s2ar1', '6s2ar2', '7ar1', '7ar2']
patchList = p['patchList']

# patchTemp = '9'
# Loop over patches 

# for i in xrange(0,2048):
# for i in xrange(0,2):


dataMapDir = '../data/'
print 'got here'

import statsTools
powers = statsTools.threedl(numCMBsets, iMax - iMin, len(patchList))

binnedPower = statsTools.onedl(len(patchList))

import pickle


for iii in xrange(iMin,iMax):
    

    noiseDirRoot = p['noiseDirRoot']


    for cmbSet in numpy.arange(0, numCMBsets):
        noiseDir = '%s_set%02d_%05d'%(noiseDirRoot, cmbSet, iii)
        try:
            os.mkdir(noiseDir)
        except:
            pass


    
        for patchNum, patch in enumerate(p['patchList']):
            print 'rank ', rank, 'of ', size , ', iii = ', iii, ' within %d to %d, patch =', patch
            tNoise = pickle.load(open(noisePsdDir+'noisePowerIAlt_'+patch+'.pkl'))
            qNoise = pickle.load(open(noisePsdDir+'noisePowerQAlt_'+patch+'.pkl'))
            uNoise = pickle.load(open(noisePsdDir+'noisePowerUAlt_'+patch+'.pkl'))


            maskM = True
            mask = liteMap.liteMapFromFits(noisePsdDir + 'weightMap_' + patchList[patchNum] + '.fits')

            # mask = liteMap.liteMapFromFits(weightMapDir + 'weightMap_'+patch+'.fits')
            # print 'b', mask.info()
            loc = numpy.where(mask.data == 0.)
            loc2 = numpy.where(mask.data < 0.)
            mask.data = numpy.sqrt(mask.data)

            # T = liteMap.liteMapFromFits(weightMapDir + 'weightMap_' + patch + '.fits')



            # if i == 0:
            #     print 'cutmap', T.info()
            TF = mask.copy()
            QF = mask.copy()
            UF = mask.copy()
            TF.data *= 0.
            QF.data *= 0.
            UF.data *= 0.

            TF.fillWithGRFFromTemplate(tNoise,bufferFactor=1)
            QF.fillWithGRFFromTemplate(qNoise,bufferFactor=1)
            UF.fillWithGRFFromTemplate(uNoise,bufferFactor=1)

            if maskM:
                TF.data /= mask.data
                QF.data /= mask.data
                UF.data /= mask.data
                TF.data[loc] = 0.
                QF.data[loc] = 0.
                UF.data[loc] = 0.
                TF.data[loc2] = 0.
                QF.data[loc2] = 0.
                UF.data[loc2] = 0.

            TF.writeFits(noiseDir+'/noise_2dPSD_v2_T_%s.fits'%(patch), overWrite=True)
            QF.writeFits(noiseDir+'/noise_2dPSD_v2_Q_%s.fits'%(patch), overWrite=True)
            UF.writeFits(noiseDir+'/noise_2dPSD_v2_U_%s.fits'%(patch), overWrite=True)

            # lower, upper, center, binnedPower[patchNum], bin_stds, bincount = statsTools.aveBinInAnnuli(tNoise)
            
            # powers[cmbSet][ iii][ patchNum] = statsTools.quickPower(TF)

            
# pickle.dump(powers, open('../data/noisePowers.pkl', 'w'))


        #for i in xrange(0,1):

                # T = liteMap.liteMapFromFits('/scratch2/r/rbond/engelen/crosssims/data/lensedCMBMaps_%05d/order5_lensedCMB_T_beam_cutout_'%(i)+ patch + '.fits')
                # if addForegrounds:
                #     FG = liteMap.liteMapFromFits('/scratch2/r/rbond/engelen/crosssims/data/noiseMaps_%03d/fg_T_beam_cutout_'%(i)+patch+'.fits')
                #     T.data += FG.data
                # Q = liteMap.liteMapFromFits('/scratch2/r/rbond/engelen/crosssims/data/lensedCMBMaps_%05d/order5_lensedCMB_Q_beam_cutout_'%(i) + patch + '.fits')
                # U = liteMap.liteMapFromFits('/scratch2/r/rbond/engelen/crosssims/data/lensedCMBMaps_%05d/order5_lensedCMB_U_beam_cutout_'%(i) + patch + '.fits')

            
