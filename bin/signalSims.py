from __future__ import print_function
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from actsims.util import seed_tracker as seedgen

from pixell import enmap, utils , lensing
from pixell import powspec, curvedsky
import numpy as np
import healpy
import subprocess, sys
from mpi4py import MPI
import time
import yaml
from enlib import bench
from pixell import aberration


def mpiMinMax(comm, iStop, iStart = 0):
#copied/pasted from Alex's "aveTools" library - https://github.com/ajvanengelen/aveTools
    rank = comm.Get_rank()
    size = comm.Get_size()

    delta = (iStop - iStart)/size
    if delta == 0:
        raise ValueError('Too many processors for too small a  loop!')

    iMin = iStart+rank*delta
    iMax = iStart+(rank+1)*delta

    if iMax>iStop:
        iMax = iStop
    elif (iMax > (iStop - delta)) and iMax <iStop:
        iMax = iStop

    return iMin, iMax, delta, rank, size


def usefulInfo(scriptFilename, inputFilename, parameterDict, gitVersion = None, runTime = None, size = None):

    output = '\n'
    output += 'This run of %s was invoked with config file %s\n\n' % (scriptFilename, inputFilename)

    if gitVersion is not None:
        output += "Git version number: %s  \n\n" % gitVersion
    output += "Parameters were:\n\n"
    for key in parameterDict.keys():
        output += "%s = %s \n" % (key, str(parameterDict[key]))
    if runTime is not None:
        output += '\nThe run took %f seconds.\n' % runTime
    if size is not None:
        output += '\nThe run was done on %i MPI processes.' % size
    return output

if __name__ == '__main__':

    try:
        with open('../inputParams/paths_local.yml') as f:
            p = yaml.load(f)
    except:
        print("ERROR: ../inputParams/paths_local.yml not found. Please copy ../inputParams/paths.yml to this file and edit with your local paths.")
        sys.exit(1)

    cmbDir = p['dataDir']

    with open('../inputParams/' + sys.argv[1]) as f:
        p = yaml.load(f)


    startTime = time.clock()




    iMin, iMax, delta, rank, size = mpiMinMax(MPI.COMM_WORLD, p['iStop'])

    shape, wcs = enmap.fullsky_geometry(p['PIX_SIZE']*utils.arcmin)
    ps = powspec.read_camb_full_lens(p['inputSpecRoot'] + "_lenspotentialCls.dat")
    lPs = powspec.read_spectrum(p['inputSpecRoot'] + "_lensedCls.dat")

    with bench.show("init"):
        Ab = aberration.Aberrator(shape, wcs, modulation=None)


    #make phi totally uncorrelated with both T and E.  This is necessary due to the way that separate phi and CMB seeds were put forward in an update to the pixell library around mid-Nov 2018
    ps[0, 1:, :] = 0.
    ps[1:, 0, :] = 0.


    start = time.time()

    for cmbSet in range(p['START_CMB_SET'], p['STOP_CMB_SET']):    
        for phiSet in range(p['START_PHI_SET'], p['STOP_PHI_SET']):    
            for iii in range(int(iMin), int(iMax)):
                print('rank', rank, 'doing cmbSet', cmbSet, 'iii' , iii, \
                    ', iMin', iMin, ', iMax', iMax, 'calling lensing.rand_map', time.time() - start)


                phiSeed = seedgen.get_phi_seed(phiSet, iii)
                cmbSeed = seedgen.get_cmb_seed(cmbSet, iii)

                with bench.show("lensing"):
                    lTquMap, = lensing.rand_map((3,)+shape, wcs, ps,
                                                              lmax=p['LMAX'],
                                                              output="l",
                                                              verbose=True,
                                                              phi_seed = phiSeed,
                                                              seed = cmbSeed)



                mapList = [lTquMap]

                mapNameList = ['fullskyLensedUnabberatedCMB']
                #Yes, there is a spelling mistake in "Unabberated" - we are
                #keeping it for reverse compatibility with previous versions
                #:(

                if p['doAberration']:


                    print('doing aberration')
                    print('rank', rank, 'doing cmbSet', cmbSet, 'iii' , iii, \
                        ', iMin', iMin, ', iMax', iMax, 'calling aberration.boost_map', time.time() - start)


                    #Note - we are aberrating and not modulating! The
                    #modulation is a frequency-dependent, so is done
                    #later.
                    with bench.show("boost"):
                        lTquMapAberrated = Ab.aberrate(lTquMap)

                    mapList += [lTquMapAberrated]
                    mapNameList += ['fullskyLensedAbberatedCMB']
                    #Yes, there is a spelling mistake in "Unabberated" - we are
                    #keeping it for reverse compatibility with previous versions
                    #:(


                for mi, mmm in enumerate(mapList):
                    print(iii, ' calling curvedsky.map2alm for ', mapNameList[mi])
                    alm = curvedsky.map2alm(mmm, lmax=p['LMAX_WRITE'])


                    filename = cmbDir + "/%s_alm_%s%s%05d.fits" \
                               % ( mapNameList[mi],
                                   ('cmbset%02d_' % cmbSet if 'CMB' in mapNameList[mi] else '' ) ,
                                   ('phiset%02d_' % phiSet \
                                    if (('Lensed' in mapNameList[mi]) or ('Phi' in mapNameList[mi])) \
                                    else '' ) ,
                                   iii)

                    print('writing to disk, filename is', filename)

                    healpy.fitsfunc.write_alm(filename ,
                                               np.complex64(alm), overwrite = True)


    if rank == 0:
        gitVersion = subprocess.check_output("git rev-parse --short HEAD", shell=True).rstrip("\n")
        with open(cmbDir + "/__README.txt", 'w') as f:
            print(usefulInfo(__file__, sys.argv[1], p, gitVersion, time.time() - start, size), file = f)

