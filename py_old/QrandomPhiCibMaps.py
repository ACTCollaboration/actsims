
#this was an attempt to fix something by using sudeeps installation - turned out not necessary
# import sys
# try:
#     sys.path.remove('/scratch/r/rbond/engelen/flipper_sudeep/python')
# except:
#     pass
# sys.path.append('/home/r/rbond/sudeep/python/flipper/python/')
from flipper import *
from mpi4py import MPI


from flipperPol import *


import sys
sys.path.append('/scratch/r/rbond/engelen/aveTools/') 
sys.path.append('/scratch2/r/rbond/engelen/limber/py/')
from limberTools import *
def phiToKappa(phiMap):

    kappaMap = phiMap.copy()

    phiF = fftTools.fftFromLiteMap(phiMap)
    
    kappaF = phiF.copy()
    kappaF.kMap *= kappaF.modLMap**2 / 2.

    kappaMap.data[:] = kappaF.mapFromFFT( setMeanToZero = True)



    return kappaMap

doFigs = False



p = flipperDict.flipperDict()
p.readFromFile(sys.argv[1])


iStart = p['iStart']
iStop = p['iStop']
phiDirRoot = p['phiDirRoot']

# The different deep fields
mapFiles= p['mapFiles']

#The position on the sky of the fields.
Ra0Array= p['Ra0Array']
Ra1Array= p['Ra1Array']
Dec0Array = p['Dec0Array']
Dec1Array = p['Dec1Array']

#buffer: buffer !=0 make the map non periodic
buffer= p['buffer']

seedbase1 = p['seedbase1']
nsNames = p['nsNames']
# MPI implementation

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

print(rank, size)

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


#Read the phi power spectrum
X = numpy.loadtxt(p['theoryScal'])
tcmb = 2.726e6
lphi = X[:,0]
clphi = X[:,4]/(lphi**4)/tcmb**2

(l_limber, clphi_limber, clcib_limber, clcibphi_limber)  = numpy.transpose(numpy.loadtxt('../../limber/data/cib_kappa_powers_non-linear.txt'))

from scipy import interpolate

f_cib = interpolate.interp1d(l_limber, clcib_limber, bounds_error = False, fill_value = 0.)
f_cibphi = interpolate.interp1d(l_limber, clcibphi_limber , bounds_error = False, fill_value = 0.)

cibPower = f_cib(lphi) 
cibPhiPower = f_cibphi(lphi) 

if doFigs:
    figure(3)
    plot(lphi,clphi * lphi**4, 'r')
    plot(l_limber,clphi_limber * l_limber**4, 'b')
   
    show()

    
globalnum = 10 * iMin + seedbase1
# for iii in xrange(iMin,1):

for iii in xrange(iMin,iMax):
    phiDir = "%s_%05d"%(phiDirRoot,iii)
    print('iteration' ,iii, iMin, iMax)
    try:
        os.mkdir(phiDir)
    except:
        pass

    count = 0 

    for i, patch in enumerate(nsNames):
        print("reading template %03d"%count)
        m = liteMap.liteMapFromFits(mapFiles[i])
        m.data[:] = 0.


        globalnum += 1 

        # this is the way of not using a template.
        m = liteMap.makeEmptyCEATemplateAdvanced(Ra0Array[i] - buffer, Dec0Array[i] - buffer, \
                                         Ra1Array[i] + buffer, Dec1Array[i] + buffer)


        print(globalnum)
        # (phiMap, cibMap) = simCorrelatedMaps(m,l_limber,clphi_limber,clcib_limber,clcibphi_limber, seed = globalnum)        
        (phiMap, cibMap) = simCorrelatedMaps(m,lphi,clphi, cibPower,cibPhiPower, seed = globalnum)        

        kappaMap = phiToKappa(phiMap)

        if doFigs:
            figure(0)
            clf()
            imshow(kappaMap.data)
            title(r'$\kappa$')
            colorbar()

            figure(1)
            clf()
            imshow(cibMap.data)
            title(r'CIB (Jy/Sr)')
            colorbar()


        
        # select the relevant part
        # m=m.selectSubMap(Ra0Array[count]-buffer,Ra1Array[count]+buffer,Dec0Array[count]-buffer,Dec1Array[count]+buffer)
        
        print('='*10)
        print('Map Coordinates')
        phiMap.info()
        print('='*10)
        
        print('='*10)
        print('Map Coordinates')
        cibMap.info()
        print('='*10)
        
        
        # fill it with a Gaussian Random realization of Phi
        # p.fillWithGaussianRandomField(l_limber,clphi_limber,bufferFactor=1)


        # take out the mean
        cibMap.data[:] -= cibMap.data.mean()
        phiMap.data[:] -= phiMap.data.mean()

        phiMap.writeFits('%s/phiMap_%s.fits'%(phiDir, patch),overWrite=True)
        kappaMap.writeFits('%s/kappaMap_%s.fits'%(phiDir, patch),overWrite=True)
        cibMap.writeFits('%s/cibMap_%s.fits'%(phiDir, patch),overWrite=True)

        print('')
        count += 1

        #this was just to check the cutout coords.
        if False:
            cutoutphi =    phiMap.selectSubMap(Ra0Array[i],Ra1Array[i],Dec0Array[i],Dec1Array[i])
            print('='*10)
            print('Cutout Coordinates')
            cutoutphi.info()
            print('='*10)

    
