
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


import sys
sys.path.append('../../')

from pixell import enmap, utils, lensing, powspec, curvedsky
import numpy as np
import healpy
import matplotlib.pyplot as plt
import os
from mpi4py import MPI
sys.path.append('../../aveTools/')
import aveTools
import pickle
import cmblens.flipper.flipperDict

# p = flipper.flipperDict.flipperDict()
p = cmblens.flipper.flipperDict.flipperDict()

p.read_from_file('../inputParams/' + sys.argv[1])

import time
startTime = time.clock()



iMin, iMax, delta, rank, size = aveTools.mpiMinMax(MPI.COMM_WORLD, p['iStop'])


# import pdb
def tqu2teb(tqu, LMAX, wantCl = False):
    alm = curvedsky.map2alm(tqu, lmax=LMAX)
    teb = curvedsky.alm2map(alm[:,None], tqu.copy()[:,None], spin=0)[:,0]
    if wantCl:
        cls = healpy.sphtfunc.alm2cl(alm)
        return teb, cls
    else:
        return teb

def phi2kappa(phiMap, LMAX):

    phiAlm = curvedsky.map2alm(phiMap, lmax = LMAX)
    ells = np.arange(LMAX-1)
    kappaAlm = healpy.sphtfunc.almxfl(phiAlm, ells * (ells + 1) / 2.)
    kappaMap = curvedsky.alm2map(kappaAlm, phiMap.copy() )
    return kappaMap




#######################

# LMAX = 5000
# LMAX_NYQ = 5400
# PIX_SIZE = 2.0

# doAll = True
# scratch = "/global/cscratch1/sd/engelen/"
# dataDir = '%s/simsS1516/data/' % scratch

# inputSpecRoot = '../input/cosmo2017' #erminia cosmology



shape, wcs = enmap.fullsky_geometry(p['PIX_SIZE']*utils.arcmin)
ps = powspec.read_camb_full_lens(p['inputSpecRoot'] + "_lenspotentialCls.dat")
lPs = powspec.read_spectrum(p['inputSpecRoot'] + "_lensedCls.dat")

doAll = True    
cmbSet = 0 # still to-do: loop over sets.

if doAll:
    uTebCls = aveTools.onedl(p['iStop'] - p['iStart'])
    lTebCls = aveTools.onedl(p['iStop'] - p['iStart'])


    for iii in range(iMin, iMax):
        print 'rank', rank, 'doing iii' , iii, ', iMin', iMin, ', iMax', iMax

        # uTquMap, lTquMap, pMap, kappaMap = lensing.rand_map((3,)+shape, wcs, ps, lmax=p['LMAX'], output="ulpk", verbose=True)

        #this is the good one! testing a change.
        
        # uTquMap, lTquMap, pMap = lensing.rand_map((3,)+shape, wcs, ps, lmax=p['LMAX'], output="ulp", verbose=True)
        uTquMap, lTquMap, pMap = lensing.rand_map((3,)+shape, wcs, ps, lmax=p['LMAX'], output="ulp", verbose=True,
                                                  separate_phi_from_cmb = True,
                                                  phi_seed = iii,
                                                  seed = iii * 100)
        

        if p['doAberration']:
            unaberrated = lTquMap.copy()

            from enlib import aberration
            lTquMap = aberration.aberrate(lTquMap,
                                          aberration.dir_equ,
                                          aberration.beta, modulation = False)

        uTebMap, uTebCls[iii] = tqu2teb(uTquMap,p['LMAX_NYQ'], wantCl = True)
        lTebMap, lTebCls[iii] = tqu2teb(lTquMap,p['LMAX_NYQ'], wantCl = True)
        

        

        cmbDir = p['dataDir'] 
        # try:
        #     os.mkdir(cmbDir)
        #     os.chmod(cmbDir, 0755)
        # except:
        #     pass

        for x in range(0, 3):
            enmap.write_map(cmbDir \
                            + "cmb_set%02d_fullskyUnlensedMap_%s_%05d.fits" % (  cmbSet, 'TQU'[x], iii),
                            np.float32(uTquMap[x]))

            enmap.write_map(cmbDir \
                            + "cmb_set%02d_fullskyLensedMap_%s_%05d.fits" % (  cmbSet, 'TQU'[x], iii),
                            np.float32(lTquMap[x]))


            os.chmod(cmbDir \
                     + "cmb_set%02d_fullskyUnlensedMap_%s_%05d.fits" % (  cmbSet, 'TQU'[x], iii), 0644)

            os.chmod(cmbDir \
                     + "cmb_set%02d_fullskyLensedMap_%s_%05d.fits" % (  cmbSet, 'TQU'[x], iii), 0644)

            if p['doAberration']:
                enmap.write_map(cmbDir \
                            + "cmb_set%02d_fullskyLensedUnaberratedMap_%s_%05d.fits" % (  cmbSet, 'TQU'[x], iii),
                            np.float32(unaberrated[x]))

                os.chmod(cmbDir \
                            + "cmb_set%02d_fullskyLensedUnaberratedMap_%s_%05d.fits" % (  cmbSet, 'TQU'[x], iii),0644)

        if False:
            for x in range(1, 3):
                enmap.write_map(cmbDir + "fullskyUnlensedMap_%s_%05d.fits" % ( 'TEB'[x], iii), np.float32(uTebMap[x]))
                enmap.write_map(cmbDir + "fullskyLensedMap_%s_%05d.fits" % ( 'TEB'[x], iii), np.float32(lTebMap[x]))




        # phiDir = p['dataDir'] + 'phi_%05d/'% (  iii)
        # try:
        #     os.mkdir(phiDir)
        # except:
        #     pass


        if False:
            enmap.write_map(cmbDir + "phiMap_%05d.fits" % iii,  pMap)
        if True:
            kappaMap = phi2kappa(pMap, LMAX = p['LMAX'])

            enmap.write_map(cmbDir + "kappaMap_%05d.fits" % iii,  np.float32(kappaMap))

            # enmap.write_map(phiDir + "kappaMap_%05d.fits" % iii,  kappaMap)

    aveTools.mpiSendReceiveList(uTebCls, MPI.COMM_WORLD, iMin, iMax, delta)
    aveTools.mpiSendReceiveList(lTebCls, MPI.COMM_WORLD, iMin, iMax, delta)


    if rank == 0:
        pickle.dump(uTebCls, open(p['dataDir'] +  'uTebClsFullsky.pkl', "wb"))
        pickle.dump(lTebCls, open(p['dataDir'] +  'lTebClsFullsky.pkl', "wb"))
    print 'finished, took' , time.clock() - startTime



    exit()
# else:
#     uTebCls = pickle.dump(uTebCls, open(dataDir + 'uTebCls.dat', "wb"))
#     pickle.dump(lTebCls, open(dataDir + 'uTebCls.dat', "wb"))
#     print 'finished, took' , time.clock() - startTime



                                            # pickle.dump(powers, open(p['workDir'] + p['basename'] + 'PowersSandbox.pkl', "wb"))

    
    stop                                                
