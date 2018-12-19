

#!/usr/bin/env python
from flipper import *
from flipperPol import *
from mpi4py import MPI
import scipy

import statsTools





def lensMaps(phi,T_map,Q_map,U_map):#,iii,count): van Engelen commented out the last two args as they don't seem to be used in this routine.
    
    
    def binomial(n,k):
        "Compute n factorial by a direct multiplicative method"
        if k > n-k: k = n-k  # Use symmetry of Pascal's triangle
        accum = 1
        for i in range(1,k+1):
            accum *= (n - (k - i))
            accum /= i
        return accum

    ft=fftTools.fftFromLiteMap(phi)
    lx_array = numpy.zeros([ft.Ny,ft.Nx])
    ly_array = numpy.zeros([ft.Ny,ft.Nx])
    
    for i in range(ft.Ny):
        lx_array[i,:]=ft.lx
    for i in range(ft.Nx):
        ly_array[:,i]=ft.ly
    
    alphaX=numpy.real(numpy.fft.ifft2(1j*lx_array*ft.kMap))
    alphaY=numpy.real(numpy.fft.ifft2(1j*ly_array*ft.kMap))

    iy,ix = numpy.mgrid[0:phi.Ny,0:phi.Nx]
    
    alphaX0 = numpy.array(numpy.round(alphaX/ T_map.pixScaleX),dtype='int64')
    alphaY0 = numpy.array(numpy.round(alphaY/ T_map.pixScaleY),dtype='int64')

    delta_alphaX=alphaX-alphaX0*T_map.pixScaleX
    delta_alphaY=alphaY-alphaY0*T_map.pixScaleY

    lensed_T_Map = T_map.copy()
    lensed_Q_Map = Q_map.copy()
    lensed_U_Map = U_map.copy()
    cont= T_map.copy()

    lensed_T_Map.data[:]=T_map.data[(iy+alphaY0)%T_map.Ny, (ix+alphaX0)%T_map.Nx]
    lensed_Q_Map.data[:]=Q_map.data[(iy+alphaY0)%T_map.Ny, (ix+alphaX0)%T_map.Nx]
    lensed_U_Map.data[:]=U_map.data[(iy+alphaY0)%T_map.Ny, (ix+alphaX0)%T_map.Nx]

    ft=fftTools.fftFromLiteMap(T_map)
    fQ=fftTools.fftFromLiteMap(Q_map)
    fU=fftTools.fftFromLiteMap(U_map)

    for n in range(1,TaylOrder):
        
        cont.data[:]=0
        for k in range(n+1):
                
            print(k, n-k, binomial(n,k), scipy.misc.factorial(n))
                
            fac=1j**n*binomial(n,k)*lx_array**(n-k)*ly_array**k/(scipy.misc.factorial(n))
            T_add=numpy.real(numpy.fft.ifft2(fac*ft.kMap))[(iy+alphaY0)%T_map.Ny, (ix+alphaX0)%T_map.Nx]*delta_alphaX**(n-k)*delta_alphaY**k
            Q_add=numpy.real(numpy.fft.ifft2(fac*fQ.kMap))[(iy+alphaY0)%T_map.Ny, (ix+alphaX0)%T_map.Nx]*delta_alphaX**(n-k)*delta_alphaY**k
            U_add=numpy.real(numpy.fft.ifft2(fac*fU.kMap))[(iy+alphaY0)%T_map.Ny, (ix+alphaX0)%T_map.Nx]*delta_alphaX**(n-k)*delta_alphaY**k
                
            lensed_T_Map.data[:]+=T_add
            lensed_Q_Map.data[:]+=Q_add
            lensed_U_Map.data[:]+=U_add
            
            cont.data[:]+=T_add



    return(lensed_T_Map, lensed_Q_Map,lensed_U_Map)



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


Ra0Array= p['Ra0Array']
Ra1Array=  p['Ra1Array']
Dec0Array =  p['Dec0Array']
Dec1Array =  p['Dec1Array']
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
epsilon = p['epsilon']

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

    for cmbSet in numpy.arange(0,numCMBsets):
        lensDir = '%s_set%02d_%05d'%(lensDirRoot, cmbSet, iii)
        unlensDir = '%s_set%02d_%05d'%(unlensDirRoot, cmbSet, iii)
        noiseDir = '%s_set%02d_%05d'%(noiseDirRoot, cmbSet, iii)


        try:
            os.mkdir(lensDir)
            os.mkdir(unlensDir)
            os.mkdir(noiseDir)
        except:
            pass

        
        # for count in range(nPatches):
        for i, patch in enumerate(nsNames):

            
            

            phi = liteMap.liteMapFromFits(phiDir+os.path.sep+'phiMap_%s.fits'%(patch))

            phiCutout = phi.selectSubMap(Ra0Array[i] - epsilon,\
                                             Ra1Array[i] + epsilon, \
                                             Dec0Array[i] + epsilon,\
                                             Dec1Array[i] - epsilon)
            stop
            T_map,Q_map,U_map =liteMapPol.simPolMapsFromEandB(phi,l,cl_TT,cl_EE,cl_TE,cl_BB,fullBeamMatrix=fullBeamMatrix,beam1d=beam1d_off)

            T_map_cutout = T_map.selectSubMap(Ra0Array[i],Ra1Array[i],Dec0Array[i],Dec1Array[i])

            stop
