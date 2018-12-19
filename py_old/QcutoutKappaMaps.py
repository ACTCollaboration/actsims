#!/usr/bin/env python
from flipper import *
from flipperPol import *
from mpi4py import MPI
import scipy

import statsTools







print("Reading dict file")
p = flipperDict.flipperDict()
p.read_from_file(sys.argv[1])


#MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
    
print(rank, size)

iStart = p['iStart']
iStop = p['iStop']
    
delta = (iStop - iStart)/size

if delta == 0:
	raise ValueError, 'Too many processors for too small a  loop!'

print(delta)
iMin = iStart+rank*delta
iMax = iStart+(rank+1)*delta

if iMax>iStop:
	iMax = iStop
elif (iMax > (iStop - delta)) and iMax <iStop:
	iMax = iStop


#Read the data power spectrum
theoryPower = numpy.loadtxt(p['theoryScal'])
l=theoryPower[:,0]
cl_TT=theoryPower[:,1]
cl_EE=theoryPower[:,2]
cl_TE=theoryPower[:,3]
cl_BB=None

#for unlensed maps with lensed power
theoryPower_lensed = numpy.loadtxt(p['theoryLens'])
l_len=theoryPower_lensed[:,0]
cl_TT_len=theoryPower_lensed[:,1]
cl_EE_len=theoryPower_lensed[:,2]
cl_TE_len=theoryPower_lensed[:,3]
cl_BB=None

lMax = 9000 if 9000 < max(l_len) else max(l_len)

l_len=l_len[:lMax]
cl_TT_len=cl_TT_len[:lMax]*2*numpy.pi/(l_len*(l_len+1))
cl_EE_len=cl_EE_len[:lMax]*2*numpy.pi/(l_len*(l_len+1))
cl_TE_len=cl_TE_len[:lMax]*2*numpy.pi/(l_len*(l_len+1))

nsNames = p['nsNames']
nsNamesSubarr = p['nsNamesSubarr']

#noise power

noiseLevels = p['noiseLevels']

cl_TT_noises = statsTools.twodl(len(noiseLevels), len(noiseLevels[0]))
cl_EE_noises = statsTools.twodl(len(noiseLevels), len(noiseLevels[0]))
cl_BB_noises = statsTools.twodl(len(noiseLevels), len(noiseLevels[0]))
cl_TE_noises = statsTools.twodl(len(noiseLevels), len(noiseLevels[0]))


for i, patch in enumerate(nsNames):
    for n, noise in enumerate(noiseLevels[i]):
        noise_ster = (numpy.pi / (180. * 60))**2 * noise**2



        cl_TT_noises[i][n] = numpy.empty(lMax)
        cl_TT_noises[i][n].fill(noise_ster)

        cl_EE_noises[i][n] = numpy.empty(lMax)
        cl_EE_noises[i][n].fill(noise_ster * 2.)

        cl_BB_noises[i][n] = numpy.empty(lMax)
        cl_BB_noises[i][n].fill(noise_ster * 2.)

        cl_TE_noises[i][n] = numpy.empty(lMax)
        cl_TE_noises[i][n].fill(0.)


(ls_fg, foregroundPower) = numpy.loadtxt('../../limber/data/foreground_powers.txt').transpose()






l=l[:lMax]
cl_TT=cl_TT[:lMax]*2*numpy.pi/(l*(l+1))
cl_EE=cl_EE[:lMax]*2*numpy.pi/(l*(l+1))
cl_TE=cl_TE[:lMax]*2*numpy.pi/(l*(l+1))


Ra0Array = p['Ra0Array']
Ra1Array =  p['Ra1Array']
Dec0Array =  p['Dec0Array']
Dec1Array =  p['Dec1Array']


# epsilonRa0Array = p['epsilonRa0Array']
# epsilonRa1Array = p['epsilonRa1Array']
# epsilonDec0Array =  p['epsilonDec0Array']
# epsilonDec1Array =  p['epsilonDec1Array']


# for i, patch in enumerate(nsNames):

#     Ra0Array[i] += epsilonRa0Array[i]
#     Ra1Array[i] += epsilonRa1Array[i]
#     Dec0Array[i] += epsilonDec0Array[i]
#     Dec1Array[i] += epsilonDec1Array[i]



Ra0ArraySubarr = p['Ra0ArraySubarr']
Ra1ArraySubarr =  p['Ra1ArraySubarr']
Dec0ArraySubarr =  p['Dec0ArraySubarr']
Dec1ArraySubarr =  p['Dec1ArraySubarr']

epsilonRa0ArraySubarr = p['epsilonRa0ArraySubarr']
epsilonRa1ArraySubarr = p['epsilonRa1ArraySubarr']
epsilonDec0ArraySubarr =  p['epsilonDec0ArraySubarr']
epsilonDec1ArraySubarr =  p['epsilonDec1ArraySubarr']



buffer=p['buffer']

phiDirRoot = p['phiDirRoot']
unlensDirRoot = p['unlensDirRoot']
lensDirRoot = p['lensDirRoot']
noiseDirRoot = p['noiseDirRoot']


templates=p['mapFiles']
TaylOrder= p['TaylOrder']

seedbase2 = p['seedbase2']

beam1d_off = {'apply':False}

# beam1d_on = {'apply':True, 'file':'beams_10000_AR2_2010_120422.dat'}


numCMBsets = p['numCMBsets']


for i, patch in enumerate(nsNames):
    for n, subarrName in enumerate(nsNamesSubarr[i]):

        Ra0ArraySubarr[i][n] += epsilonRa0ArraySubarr[i][n]
        Ra1ArraySubarr[i][n] += epsilonRa1ArraySubarr[i][n]
        Dec0ArraySubarr[i][n] += epsilonDec0ArraySubarr[i][n]
        Dec1ArraySubarr[i][n] += epsilonDec1ArraySubarr[i][n]


            



fullBeamMatrix = {'apply':False}


# beamDir = '/scratch2/r/rbond/engelen/actmaps/data/beams_140206/'
# beamName = 'beam_tform_140206_jitmap_deep6.txt'
# beamFile = beamDir + beamName

# beamData = (numpy.loadtxt(beamFile)).transpose()


nPatches=len(templates)

globalnum = 10 * iMin + seedbase2

# dataMapDir = '/scratch2/r/rbond/engelen/lensRecon/maps/dataMaps/actpolDeep/'
#changing for local beams files...
dataMapDir = '../data/'
print('got here')



for iii in xrange(iMin,iMax):
#for iii in xrange(iMin,1):
    globalnum += 1

    phiDir = '%s_%05d'%(phiDirRoot,iii)

    for i, patch in enumerate(nsNames):

        phi = liteMap.liteMapFromFits(phiDir+os.path.sep+'phiMap_%s.fits'%(patch))            
        kappa = liteMap.liteMapFromFits(phiDir+os.path.sep+'kappaMap_%s.fits'%(patch))            

        for n, subarrName in enumerate(nsNamesSubarr[i]):

            print('*** iter ', iii, 'from ', iMin, 'to ', iMax , \
                ' -- patch ', patch, i, ' of ', len(nsNames), \
                ' -- subarr ' , subarrName, n, ' of ', len(nsNamesSubarr[i]))

            kappa_cutout = kappa.selectSubMap(Ra0ArraySubarr[i][n],Ra1ArraySubarr[i][n],
                                              Dec0ArraySubarr[i][n],Dec1ArraySubarr[i][n])

            
            kappa_cutout.writeFits(phiDir + os.path.sep + 'kappaMap_cutout_%s_%s.fits' % (patch, subarrName), overWrite = True)


stop            

    # for cmbSet in numpy.arange(0,numCMBsets):
    #     lensDir = '%s_set%02d_%05d'%(lensDirRoot, cmbSet, iii)
    #     unlensDir = '%s_set%02d_%05d'%(unlensDirRoot, cmbSet, iii)
    #     noiseDir = '%s_set%02d_%05d'%(noiseDirRoot, cmbSet, iii)


        # try:
        #     os.mkdir(lensDir)
        #     os.mkdir(unlensDir)
        #     os.mkdir(noiseDir)
        # except:
        #     pass

        
        # for count in range(nPatches):

            
