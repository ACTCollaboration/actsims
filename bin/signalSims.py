
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from actsims.util import seed_tracker as seedgen

import sys
sys.path.append('../../')

from pixell import enmap, utils , lensing
from pixell import powspec, curvedsky
import numpy as np
import healpy
import matplotlib.pyplot as plt
import os
from mpi4py import MPI
import pickle, time
from actsims import flipperDict



def mpiMinMax(comm, iStop, iStart = 0):
#copied/pasted from Alex's "aveTools" library - https://github.com/ajvanengelen/aveTools
    rank = comm.Get_rank()
    size = comm.Get_size()

    delta = (iStop - iStart)/size
    if delta == 0:
        raise ValueError, 'Too many processors for too small a  loop!'

    iMin = iStart+rank*delta
    iMax = iStart+(rank+1)*delta

    if iMax>iStop:
        iMax = iStop
    elif (iMax > (iStop - delta)) and iMax <iStop:
        iMax = iStop

    return iMin, iMax, delta, rank, size



# p = flipper.flipperDict.flipperDict()
p = flipperDict.flipperDict()

p.read_from_file('../inputParams/' + sys.argv[1])

startTime = time.clock()



iMin, iMax, delta, rank, size = mpiMinMax(MPI.COMM_WORLD, p['iStop'])

shape, wcs = enmap.fullsky_geometry(p['PIX_SIZE']*utils.arcmin)
ps = powspec.read_camb_full_lens(p['inputSpecRoot'] + "_lenspotentialCls.dat")
lPs = powspec.read_spectrum(p['inputSpecRoot'] + "_lensedCls.dat")

doAll = True    

#make phi totally uncorrelated with both T and E.  This is necessary due to the way that separate phi and CMB seeds were put forward in an update to the pixell library around mid-Nov 2018
ps[0, 1:, :] = 0.
ps[1:, 0, :] = 0.


start = time.time()

for cmbSet in range(p['START_SET'], p['STOP_SET']):    

    for iii in range(iMin, iMax):
        print('rank', rank, 'doing cmbSet', cmbSet, 'iii' , iii, \
            ', iMin', iMin, ', iMax', iMax, 'calling lensing.rand_map', time.time() - start)


        phiSeed = seedgen.get_phi_seed(iii)
        cmbSeed = seedgen.get_cmb_seed(cmbSet, iii)

        uTquMap, lTquMap, pMap = lensing.rand_map((3,)+shape, wcs, ps, lmax=p['LMAX'], output="ulp", verbose=True,
                                                          phi_seed = phiSeed,
                                                          seed = cmbSeed)


                
        mapList = [uTquMap, lTquMap, pMap]

        mapNameList = ['fullskyUnlensedCMB', 'fullskyLensedUnabberatedCMB', 'fullskyPhi']

        if p['doAberration']:


            from pixell import aberration
            print('doing aberration')
            #This was Sigurd's old version
            # lTquMap = aberration.aberrate(lTquMap,
            #                               aberration.dir_equ,
            #                               aberration.beta, modulation = False)
            print('rank', rank, 'doing cmbSet', cmbSet, 'iii' , iii, \
                ', iMin', iMin, ', iMax', iMax, 'calling aberration.boost_map', time.time() - start)

            lTquMapAberrated, AForMod = aberration.boost_map(lTquMap,
                                                    aberration.dir_equ,
                                                    aberration.beta,
                                                    modulation = None,
                                                    return_modulation = True)

            mapList += [lTquMapAberrated]
            mapNameList += ['fullskyLensedCMB']
            

        for mi, mmm in enumerate(mapList):
            print(iii, ' calling curvedsky.map2alm for ', mapNameList[mi])
            alm = curvedsky.map2alm(mmm, lmax=p['LMAX_WRITE'])

            cmbDir = p['dataDir']
            print('writing to disk')


            filename = cmbDir + "/%s_alm_%s%05d.fits" \
                       % ( mapNameList[mi],
                           ('set%02d_' % cmbSet if 'CMB' in mapNameList[mi] else '' ) ,
                           iii)

            healpy.fitsfunc.write_alm(filename ,
                                       np.complex64(alm), overwrite = True)










































# def tqu2teb(tqu, LMAX, wantCl = False, wantAlmAndCl = False):
#     alm = curvedsky.map2alm(tqu, lmax=LMAX)
#     teb = curvedsky.alm2map(alm[:,None], tqu.copy()[:,None], spin=0)[:,0]
#     if wantCl:
#         cls = healpy.sphtfunc.alm2cl(alm)
#         return teb, cls
#     if wantAlmAndCl:
#         cls = healpy.sphtfunc.alm2cl(alm)
#         return teb, alm, cls
#     else:
#         return teb

# def phi2kappa(phiMap, LMAX):

#     phiAlm = curvedsky.map2alm(phiMap, lmax = LMAX)
#     ells = np.arange(LMAX-1)
#     kappaAlm = healpy.sphtfunc.almxfl(phiAlm, ells * (ells + 1) / 2.)
#     kappaMap = curvedsky.alm2map(kappaAlm, phiMap.copy() )
#     return kappaMap









            # cls[iii] += [healpy.sphtfunc.alm2cl(alm)]
                                      



                                      
        # uTebMap, uTebAlms, uTebCls[iii] = tqu2teb(uTquMap,p['LMAX_NYQ'], wantAlmAndCl = True)


        # lTebMap, uTebAlms, lTebCls[iii] = tqu2teb(lTquMap,p['LMAX_NYQ'], wantAlmAndCl = True)
        

        

        # try:
        #     os.mkdir(cmbDir)
        #     os.chmod(cmbDir, 0755)
        # except:
        #     pass
        # alm = curvedsky.map2alm(tqu, lmax=LMAX, alm2map)


        # stop
        # for x in range(0, 3):
        #     enmap.write_map(cmbDir \
        #                     + "cmb_set%02d_fullskyUnlensedMap_%s_%05d.fits" % (  cmbSet, 'TQU'[x], iii),
        #                     np.float32(uTquMap[x]))

        #     enmap.write_map(cmbDir \
        #                     + "cmb_set%02d_fullskyLensedMap_%s_%05d.fits" % (  cmbSet, 'TQU'[x], iii),
        #                     np.float32(lTquMap[x]))


        #     os.chmod(cmbDir \
        #              + "cmb_set%02d_fullskyUnlensedMap_%s_%05d.fits" % (  cmbSet, 'TQU'[x], iii), 0644)

        #     os.chmod(cmbDir \
        #              + "cmb_set%02d_fullskyLensedMap_%s_%05d.fits" % (  cmbSet, 'TQU'[x], iii), 0644)

        #     if p['doAberration']:
        #         enmap.write_map(cmbDir \
        #                     + "cmb_set%02d_fullskyLensedUnaberratedMap_%s_%05d.fits" % (  cmbSet, 'TQU'[x], iii),
        #                     np.float32(unaberrated[x]))

        #         os.chmod(cmbDir \
        #                     + "cmb_set%02d_fullskyLensedUnaberratedMap_%s_%05d.fits" % (  cmbSet, 'TQU'[x], iii),0644)

        # if False:
        #     for x in range(1, 3):
        #         enmap.write_map(cmbDir + "fullskyUnlensedMap_%s_%05d.fits" % ( 'TEB'[x], iii), np.float32(uTebMap[x]))
        #         enmap.write_map(cmbDir + "fullskyLensedMap_%s_%05d.fits" % ( 'TEB'[x], iii), np.float32(lTebMap[x]))




        # # phiDir = p['dataDir'] + 'phi_%05d/'% (  iii)
        # # try:
        # #     os.mkdir(phiDir)
        # # except:
        # #     pass


        # if False:
        #     enmap.write_map(cmbDir + "phiMap_%05d.fits" % iii,  pMap)
        # if True:
        #     kappaMap = phi2kappa(pMap, LMAX = p['LMAX'])

        #     enmap.write_map(cmbDir + "kappaMap_%05d.fits" % iii,  np.float32(kappaMap))

        #     # enmap.write_map(phiDir + "kappaMap_%05d.fits" % iii,  kappaMap)

    # aveTools.mpiSendReceiveList(uTebCls, MPI.COMM_WORLD, iMin, iMax, delta)
    # aveTools.mpiSendReceiveList(lTebCls, MPI.COMM_WORLD, iMin, iMax, delta)


    # if rank == 0:
    #     pickle.dump(uTebCls, open(p['dataDir'] +  'uTebClsFullsky.pkl', "wb"))
    #     pickle.dump(lTebCls, open(p['dataDir'] +  'lTebClsFullsky.pkl', "wb"))
    # print 'finished, took' , time.clock() - startTime




# else:
#     uTebCls = pickle.dump(uTebCls, open(dataDir + 'uTebCls.dat', "wb"))
#     pickle.dump(lTebCls, open(dataDir + 'uTebCls.dat', "wb"))
#     print 'finished, took' , time.clock() - startTime



                                            # pickle.dump(powers, open(p['workDir'] + p['basename'] + 'PowersSandbox.pkl', "wb"))

    

