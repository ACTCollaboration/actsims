from flipper import *
import pdb
import scipy


def lensMaps(phi,T_map, Q_map, U_map, TaylOrder = 5):#,iii,count): van Engelen commented out the last two args as they don't seem to be used in this routine.
    
    
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






def rotMaps(alphaMapRad, T_map, Q_map, U_map):
#note T map is unaffected
    

    if not(numpy.all(numpy.isreal(alphaMapRad))):
        raise ValueError('input rotation map should be real.')

    Qloc = Q_map.copy()
    Uloc = U_map.copy()



    plus = Qloc.data + 1j * Uloc.data

    minus = Qloc.data - 1j * Uloc.data


    plusRot =   numpy.exp(+2 * 1j * alphaMapRad.data) * plus
    minusRot =  numpy.exp(-2 * 1j * alphaMapRad.data) * minus


    #note, result should be real regardless (if the initial Q/U is, that is) --
    #we just don't want an imaginary part to be stored.
    Qloc.data = numpy.real(1./2.        * (plusRot + minusRot))
    Uloc.data = numpy.real(1./(2. * 1j) * (plusRot - minusRot))


    
    return (T_map.copy(), Qloc, Uloc)


def modMaps(tauMap, T_map, Q_map, U_map, isOpticalDepth = True):

    if not(numpy.all(numpy.isreal(tauMap))):
        raise ValueError('input mod map should be real.')

    Tloc = T_map.copy()
    Qloc = Q_map.copy()
    Uloc = U_map.copy()

    modFac = (numpy.exp(-1 * tauMap.data) if isOpticalDepth else modFac)


    Tloc.data = modFac * Tloc.data
    Qloc.data = modFac * Qloc.data
    Uloc.data = modFac * Uloc.data

    
    return (Tloc, Qloc, Uloc)
    
    
def sourceMap(sVecJy, dndsVecPerJyPerSr, templateMap, sMinJy = 1e-4, sMaxJy = 30e-3, freqGHz = 150.):

    output = templateMap.copy()

    dsJy = numpy.abs(numpy.append([0.], sVecJy[:-1] - sVecJy[1:]))

    numSourcesMean = dndsVecPerJyPerSr * dsJy * templateMap.area * (numpy.pi / 180)**2


    numSources = numpy.random.poisson(numSourcesMean)

    output.data[:] = 0.

    for si, s in enumerate(sVecJy):
        if (s < sMinJy) or (s > sMaxJy):
            continue
        for i in range(numSources[si]):
            output.data[numpy.random.randint(output.data.shape[0]), \
                            numpy.random.randint(output.data.shape[1])] \
                            += s / (templateMap.pixScaleX * templateMap.pixScaleY)

    output.convertToMicroKFromJyPerSr(freqGHz = freqGHz)

    return output



            
def phiToKappa(phiMap):

    kappaMap = phiMap.copy()

    phiF = fftTools.fftFromLiteMap(phiMap)
    
    kappaF = phiF.copy()
    kappaF.kMap *= kappaF.modLMap**2 / 2.

    kappaMap.data[:] = kappaF.mapFromFFT( setMeanToZero = True)



    return kappaMap


    
