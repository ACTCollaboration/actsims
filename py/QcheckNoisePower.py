

#!/usr/bin/env python
from flipper import *
from flipperPol import *
import scipy


import statsTools
import aveTools
import pickle

import pdb

print "Reading dict file"
p = flipperDict.flipperDict()
p.read_from_file(sys.argv[1])

start = time.time()

from mpi4py import MPI

#MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
    
print rank, size

iStart = p['iStart']
iStop = p['iStop']
    
delta = (iStop - iStart)/size

if delta == 0:
	raise ValueError, 'Too many processors for too small a  loop!'

print delta
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
cl_BB_len=theoryPower_lensed[:,3]
cl_TE_len=theoryPower_lensed[:,4]

lMax = 9000 if 9000 < max(l_len) else max(l_len)

l_len=l_len[:lMax]
cl_TT_len=cl_TT_len[:lMax]*2*numpy.pi/(l_len*(l_len+1))
cl_EE_len=cl_EE_len[:lMax]*2*numpy.pi/(l_len*(l_len+1))
cl_TE_len=cl_TE_len[:lMax]*2*numpy.pi/(l_len*(l_len+1))
cl_BB_len=cl_BB_len[:lMax]*2*numpy.pi/(l_len*(l_len+1))


nsNames = p['nsNames']
nsNamesSubarr = p['nsNamesSubarr']


#noise power

#note this was changed to noiseLevels
# noiseLevel = p['noiseLevel']
# noise_ster = (numpy.pi / (180. * 60))**2 * noiseLevel**2
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


        
cl_TT_noise = numpy.empty(lMax)
cl_TT_noise.fill(noise_ster)

cl_EE_noise = numpy.empty(lMax)
cl_EE_noise.fill(noise_ster * 2.)

cl_BB_noise = numpy.empty(lMax)
cl_BB_noise.fill(noise_ster * 2.)

cl_TE_noise = numpy.empty(lMax)
cl_TE_noise.fill(0.)


(ls_fg, foregroundPower) = numpy.loadtxt('../../limber/data/foreground_powers.txt').transpose()



l=l[:lMax]
cl_TT=cl_TT[:lMax]*2*numpy.pi/(l*(l+1))
cl_EE=cl_EE[:lMax]*2*numpy.pi/(l*(l+1))
cl_TE=cl_TE[:lMax]*2*numpy.pi/(l*(l+1))


Ra0Array= p['Ra0Array']
Ra1Array=  p['Ra1Array']
Dec0Array =  p['Dec0Array']
Dec1Array =  p['Dec1Array']
buffer = p['buffer']

phiDirRoot = p['phiDirRoot']
unlensDirRoot = p['unlensDirRoot']
lensDirRoot = p['lensDirRoot']
noiseDirRoot = p['noiseDirRoot']

numCMBsets = p['numCMBsets']


templates=p['mapFiles']
TaylOrder= p['TaylOrder']

seedbase2 = p['seedbase2']

beam1d_off = {'apply':False}

beam1d_on = {'apply':True, 'file':'beams_10000_AR2_2010_120422.dat'}


fullBeamMatrix = {'apply':False}

globalnum = 0

nPatches = len(nsNames)
tapers = [None] * nPatches
nSims = iMax - iMin



cosineApod = {'apply':True,'lenApod':200,'pad':0}

# patchNames = ['deep1', 'deep5', 'deep6']

templateMapLoc = '/scratch2/r/rbond/engelen/lensRecon/maps/dataMaps/actpolDeep/'
smallNum = 10
dataMapDir = '../data/'



cosineApod = {'apply':True,'lenApod':500,'pad':0}

# iStop = 10
iMax = 2

noisePsdDir = '/scratch2/r/rbond/engelen/new6/lensRecon/maps/dataMaps/actpolDeep/'



weightMapDir = noisePsdDir

# patchList = ['5s1ar1',  '6s1ar1', '6s2ar1', '6s2ar2', '7ar1', '7ar2']
# patchList = [p['nsNamesSubarr']][0]

doAll = True


if doAll:

    power_data =                statsTools.twodl(nPatches, len(nsNamesSubarr[0]) )
    power_noiseSims =          statsTools.fourdl(nPatches, len(nsNamesSubarr[0]), numCMBsets, iMax - iMin)

    

    
    for i, patch in enumerate(nsNames):
        for n, subarrName in enumerate(nsNamesSubarr[i]):
                
            print 'got here!'

            # dataMap = liteMap.liteMapFromFits('/scratch2/r/rbond/engelen/new4/maps/actpolDeep/dataCoadd_I_' + subarrName + '.fits')        
            power_data[i][n] = aveTools.allpowersFitsFile('/scratch2/r/rbond/engelen/new6/lensRecon/maps/dataMaps//actpolDeep/dataCoadd_', '_' + subarrName + '.fits', useI = True, TOnly = False)

# power_data[i][n] = statsTools.quickPower(dataMap, applySlepianTaper = True)

            for iii in xrange(iMin,iMax):

                print i, n, iii 



                for cmbSet in numpy.arange(0,numCMBsets):
                    noiseDir = '%s_set%02d_%05d'%(noiseDirRoot, cmbSet, iii)
                    power_noiseSims[i][n][cmbSet][iii] = aveTools.allpowersFitsFile(noiseDir + '/noise_2dPSD_', '_%s.fits' % subarrName , useI = False, TOnly = False)


                                                                         # '/scratch2/r/rbond/engelen/new4/maps/actpolDeep/dataCoadd_', '_' + subarrName + '.fits')



                    
figure('noise power')
clf()
lbin = power_data[0][0]['lbin']



i = 0

numPatches = len(nsNamesSubarr[i])
for n, subarrName in enumerate(nsNamesSubarr[i]):
    subplot(2,3, n + 1 )

    semilogy(lbin, power_data[i][n]['cl_TT'], label = 'data')

    for iii in xrange(iMin,iMax):
        for cmbSet in numpy.arange(0,numCMBsets):

            semilogy(lbin, power_noiseSims[i][n][cmbSet][iii]['cl_TT'], label = 'sim %i set %i' %(iii, cmbSet))

    title(subarrName)

legend()


show()


#         # globalnum += 1

#         phiDir = '%s_%05d'%(phiDirRoot,iii)

#         for cmbSet in numpy.arange(0,numCMBsets):
#             lensDir = '%s_set%02d_%05d'%(lensDirRoot, cmbSet, iii)
#             unlensDir = '%s_set%02d_%05d'%(unlensDirRoot, cmbSet, iii)
#             noiseDir = '%s_set%02d_%05d'%(noiseDirRoot, cmbSet, iii)






                


#         phi = liteMap.liteMapFromFits(phiDir+os.path.sep+'phiMap_%s.fits'%(patch))



#         power_unl[i][cmbSet][iii] = aveTools.allpowersFitsFile(lensDir + '/unlensedCMB_', '_%s.fits'%(patch), method = 'standard')

#                 power_lensed[i][cmbSet][iii] = aveTools.allpowersFitsFile(lensDir + '/order%d_periodic_lensedCMB_'%TaylOrder,'_%s.fits'%(patch))
                
                
#                 # T_map = liteMap.liteMapFromFits(lensDir+'/unlensedCMB_T_%s.fits'%(patch), overWrite=True)
#                 # Q_map = liteMap.liteMapFromFits(lensDir+'/unlensedCMB_Q_%s.fits'%(patch), overWrite=True)
#                 # U_map = liteMap.liteMapFromFits(lensDir+'/unlensedCMB_U_%s.fits'%(patch), overWrite=True)

#                 for iii in xrange(iMin,iMax):


#                     print '*** iter ', iii, 'from ', iMin, 'to ', iMax , \
#                         ' -- set ', cmbSet, ' of ', numCMBsets, \
#                         ' -- patch ', patch, i, ' of ', len(nsNames), \
#                         ' -- subarr ' , subarrName, n, ' of ', len(nsNamesSubarr[i])




#                     sampleCutoutMap = liteMap.liteMapFromFits(lensDir+'/order%d_lensedCMB_'%TaylOrder + 'T_beam%s_cutout_%s.fits'%(subarrName,patch))

#                     # tapers[i] = liteMapPol.initializeCosineWindow(sampleCutoutMap, cosineApod['lenApod'], cosineApod['pad'])  # taper weight map
#                     tapers[i] = None
#                     beamFile = dataMapDir + 'beam_' + subarrName + '.txt'
                



                    
#                     power_lensed_beam_cutout[i][n][cmbSet][iii] = aveTools.allpowersFitsFile(lensDir+'/order%d_lensedCMB_'%TaylOrder, '_beam%s_cutout_%s.fits'%(subarrName,patch), window = tapers[i])
#                     # T_map_lensed_beam_cutout = liteMap.liteMapFromFits(lensDir+'/order%d_lensedCMB_T_beam%s_cutout_%s.fits'%(TaylOrder,subarrName,patch), overWrite=True)
#                     # Q_map_lensed_beam_cutout = liteMap.liteMapFromFits(lensDir+'/order%d_lensedCMB_Q_beam%s_cutout_%s.fits'%(TaylOrder,subarrName,patch), overWrite=True)
#                     # U_map_lensed_beam_cutout = liteMap.liteMapFromFits(lensDir+'/order%d_lensedCMB_U_beam%s_cutout_%s.fits'%(TaylOrder,subarrName,patch), overWrite=True)


#                     # power_noise[i][n][cmbSet][iii] = aveTools.allpowersFitsFile(noiseDir+'/noise_', '_cutout_%s_%s.fits'%(subarrName, patch))
#                     # T_map_noise = liteMap.liteMapFromFits(noiseDir+'/noise_T_%s_%s.fits'%(subarrName, patch), overWrite=True)
#                     # Q_map_noise = liteMap.liteMapFromFits(noiseDir+'/noise_Q_%s_%s.fits'%(subarrName, patch), overWrite=True)
#                     # U_map_noise = liteMap.liteMapFromFits(noiseDir+'/noise_U_%s_%s.fits'%(subarrName, patch), overWrite=True)

#                     # power_T_fg_beam[i][n][cmbSet][iii] = aveTools.allpowersFitsFile(noiseDir + '/fg_','_beam%s_cutout_%s.fits'%(subarrName, patch), TOnly = True)
#                     # T_map_fg_beam = liteMap.liteMapFromFits(noiseDir + '/fg_T_beam%s_%s.fits'%(subarrName, patch), overWrite = True)

#                     if buffer !=0:

#                         # T_map_noise_cutout = liteMap.liteMapFromFits(noiseDir+'/noise_T_cutout_%s_%s.fits'%(subarrName, patch), overWrite=True)
#                         # Q_map_noise_cutout = liteMap.liteMapFromFits(noiseDir+'/noise_Q_cutout_%s_%s.fits'%(subarrName, patch), overWrite=True)
#                         # U_map_noise_cutout = liteMap.liteMapFromFits(noiseDir+'/noise_U_cutout_%s_%s.fits'%(subarrName, patch), overWrite=True)
#                         power_noise_cutout[i][n][cmbSet][iii] = aveTools.allpowersFitsFile(noiseDir + '/noise_', '_cutout_%s_%s.fits'%(subarrName, patch), window = tapers[i])


#                         power_T_fg_beam_cutout[i][n][cmbSet][iii] = aveTools.allpowersFitsFile(noiseDir + '/fg_','_beam%s_cutout_%s.fits'%(subarrName, patch), window = tapers[i])


#                         power_unl_beam_cutout[i][n][cmbSet][iii] = aveTools.allpowersFitsFile(lensDir + '/unlensedCMB_','_beam%s_cutout_%s.fits' %(subarrName, patch), window = tapers[i])

#                         # T_map_beam_cutout = liteMap.liteMapFromFits(lensDir + '/unlensedCMB_T_beam%s_cutout_%s.fits' %(subarrName, patch), overWrite=True)
#                         # Q_map_beam_cutout = liteMap.liteMapFromFits(lensDir + '/unlensedCMB_Q_beam%s_cutout_%s.fits' %(subarrName, patch), overWrite=True)
#                         # U_map_beam_cutout = liteMap.liteMapFromFits(lensDir + '/unlensedCMB_U_beam%s_cutout_%s.fits' %(subarrName, patch), overWrite=True)


#                     #unlensed maps with lensed power


#                     # power_unlLenspower_beam[i][n][cmbSet][iii] = aveTools.allpowersFitsFile(unlensDir+'/unlensedCMB_','_beam%s_%s.fits'%(subarrName, patch))
#                     # T_map_unl_beam = liteMap.liteMapFromFits(unlensDir+'/unlensedCMB_T_beam%s_%s.fits'%(subarrName, patch), overWrite=True)
#                     # Q_map_unl_beam = liteMap.liteMapFromFits(unlensDir+'/unlensedCMB_Q_beam%s_%s.fits'%(subarrName, patch), overWrite=True)
#                     # U_map_unl_beam = liteMap.liteMapFromFits(unlensDir+'/unlensedCMB_U_beam%s_%s.fits'%(subarrName, patch), overWrite=True)

#                     if buffer !=0:
#                         power_unlLenspower_beam_cutout[i][n][cmbSet][iii] = aveTools.allpowersFitsFile(unlensDir+'/unlensedCMB_','_beam%s_cutout_%s.fits'%(subarrName, patch), window = tapers[i])
#                         # T_map_unl_beam_cutout = liteMap.liteMapFromFits(unlensDir+'/unlensedCMB_T_beam%s_cutout_%s.fits'%(subarrName, patch), overWrite=True)
#                         # Q_map_unl_beam_cutout = liteMap.liteMapFromFits(unlensDir+'/unlensedCMB_Q_beam%s_cutout_%s.fits'%(subarrName, patch), overWrite=True)
#                         # U_map_unl_beam_cutout = liteMap.liteMapFromFits(unlensDir+'/unlensedCMB_U_beam%s_cutout_%s.fits'%(subarrName, patch), overWrite=True)

        

                        
                            


# print 'got here'


# def myplot(inDict, name = ''):
#     semilogy(inDict['lbin'], inDict['dlbin'], label = name, lw = 2.5, marker = 'o', ls = 'None')





# sym = '-'

# lmaxlen = len(l_len)
# n = 0
# rc('axes', color_cycle=['r', 'g', 'b', 'y', 'c', 'm', 'k', 'violet', 'orange', 'grey'])
# for i, patch in enumerate(nsNames):    
    
#     figure('allcurves periodic ' + ' '  + patch, figsize = (15,10))
#     clf()
#     myplot(power_unl[i][n][0], name = 'unl')
#     myplot(power_lensed[i][n][0], name = 'len')

#     semilogy(l,  cl_TT * l**2/2/pi, sym, label = 'unlensed')
#     # semilogy(l_len, cl_TT_len * beamData[:lmaxlen,1]**2 * l_len**2/2/pi, label = 'lensed beam')
#     semilogy(l_len, cl_TT_len * l_len**2/2/pi, label = 'lensed')

#     legend()
#     xlim([0,10e3])
#     ylim([1e-5,1e4])


#     show()
#     savefig('../plot/allcurves_periodic' + ' '  + patch + '.pdf')

#     for n, subarrName in enumerate(nsNamesSubarr[i]):
        
#         figure('allcurves ' + subarrName + ' '  + patch, figsize = (15,10))
#         clf()

#         beamFile = dataMapDir + 'beam_' + subarrName + '.txt'
#         beamData = loadtxt(beamFile)

        
#         # myplot(power_lensed_beam[i][n][0][0], 'lensed_beam')

#         myplot(power_lensed_beam_cutout[i][n][0][0], 'lensed_beam_cutout')
#         # myplot(power_noise[i][n][0][0], 'noise')
#         # myplot(power_T_fg_beam[i][n][0][0], 'T_fg_beam')
#         myplot(power_noise_cutout[i][n][0][0], 'noise_cutout')
#         myplot(power_T_fg_beam_cutout[i][n][0][0], 'T_fg_beam_cutout')
#         # myplot(power_unl_beam[0][0][0], 'unl_beam')

#         myplot(power_unl_beam_cutout[i][n][0][0], 'unl_beam_cutout')
#         myplot(power_unlLenspower_beam_cutout[i][n][0][0], 'unlLenspower_beam_cutout')
#         # myplot(power_unlLenspower_beam[i][n][0][0], 'unlLenspower_beam')

#         semilogy(l[:lmaxlen],  cl_TT[:lmaxlen] *        beamData[:lmaxlen,1]**2 * l[:lmaxlen]**2/2/pi, sym, label = 'unlensed beam')
#     # semilogy(l_len, cl_TT_len * beamData[:lmaxlen,1]**2 * l_len**2/2/pi, label = 'lensed beam')
#         semilogy(l_len, cl_TT_len * beamData[:lmaxlen,1]**2 * l_len**2/2/pi, sym, label = 'lensed beam')

        
        

#         legend()
#         xlim([0,10e3])
#         ylim([1e-5,1e4])

#         show()




#         savefig('../plot/allcurves ' + subarrName + ' '  + patch + '.pdf')

        



#                     # power_lensed_beam[i][n][cmbSet][iii] = aveTools.allpowersFitsFile(lensDir+'/order%d_periodic_lensedCMB_'%TaylOrder , \
#                     #                                                        '_beam%s_%s.fits'%(subarrName,patch))

#                     # T_map_lensed_beam = liteMap.liteMapFromFits(lensDir+'/order%d_periodic_lensedCMB_T_beam%s_%s.fits'%(TaylOrder,subarrName,patch), overWrite=True)
#                     # Q_map_lensed_beam = liteMap.liteMapFromFits(lensDir+'/order%d_periodic_lensedCMB_Q_beam%s_%s.fits'%(TaylOrder,subarrName,patch), overWrite=True)
#                     # U_map_lensed_beam = liteMap.liteMapFromFits(lensDir+'/order%d_periodic_lensedCMB_U_beam%s_%s.fits'%(TaylOrder,subarrName,patch), overWrite=True)










# # stop



# # if doAll:
# #     for i, patch in enumerate(nsNames):

# #         templateMapName = 'dataCoadd_I_' + patch + '.fits'

# #         templateMap = liteMap.liteMapFromFits(templateMapLoc + templateMapName)

# #         tapers[i] = liteMapPol.initializeCosineWindow(templateMap, cosineApod['lenApod'], cosineApod['pad'])  # taper weight map



# #     # for iii in xrange(iMin,iMax):

# #     for iii in xrange(iMin,smallNum):
# #         globalnum += 1

# #         phiDir = '%s_%05d'%(phiDirRoot,iii)

# #         lensDir = '%s_%05d'%(lensDirRoot,iii)
# #         unlensDir = '%s_%05d'%(unlensDirRoot,iii)
# #         noiseDir = '%s_%05d'%(noiseDirRoot,iii)


# #         # for count in range(nPatches):
# #         for i, patch in enumerate(nsNames):

# #             print patch, nPatches, iii, iMax - iMin, time.time() - start

# #             if buffer == 0:
# #                 print 'not implemented, just assuming for now that buffer != 0'
# #                 stop




# #             # power_len_beam_cutouts[i][iii] = aveTools.allpowersFitsFile(lensDir + '/order%d_lensedCMB_'%TaylOrder, '_beam_cutout_%s.fits'%(patch), window = tapers[i], beamFile = templateMapLoc + 'beam_' + patch + '.txt')

# #             # power_unl_beam_cutouts[i][iii] = aveTools.allpowersFitsFile(unlensDir + '/unlensedCMB_', '_beam_cutout_%s.fits'%(patch), window = tapers[i])
# #             power_len_beam_cutouts[i][iii] = aveTools.allpowersFitsFile(lensDir + '/order%d_periodic_lensedCMB_'%TaylOrder, '_beam_%s.fits'%(patch), window = None, beamFile = templateMapLoc + 'beam_' + patch + '.txt')

# #             power_unl_beam_cutouts[i][iii] = aveTools.allpowersFitsFile(unlensDir + '/unlensedCMB_', '_beam_cutout_%s.fits'%(patch), window = tapers[i])


# #             pickle.dump(power_len_beam_cutouts, open('../data/power_len_beam_cutouts.pkl', 'w'))
# #     pickle.dump(power_unl_beam_cutouts, open('../data/power_unl_beam_cutouts.pkl', 'w'))
        


# # else:

# #     power_len_beam_cutouts = pickle.load(open('../data/power_len_beam_cutouts.pkl', 'r'))
# #     power_unl_beam_cutouts = pickle.load(open('../data/power_unl_beam_cutouts.pkl', 'r'))
    

# # powers = [ 'cl_TT', 'cl_TE',  'cl_EE', 'cl_EB', 'cl_BB']
# # powersPretty = ['C_l^{TT}' ,'C_l^{TE}', 'C_l^{EE}', 'C_l^{EB}', 'C_l^{BB}']
# # theoryLensedPowers = [cl_TT_len, cl_TE_len,  cl_EE_len, zeros(len(l_len)), cl_BB_len]
# # nPowers = len(powers)
# # lbin = power_len_beam_cutouts[0][0]['lbin']

# # statsLenBeamCutouts = statsTools.twodl(nPatches, nPowers)
# # statsUnlBeamCutouts = statsTools.twodl(nPatches, nPowers)


# # for i, patch in enumerate(nsNames):

# #     for j, power in enumerate(powers):

# #         clsloc = [power_len_beam_cutouts[i][k][power] for k in arange(iMin, smallNum)]

# #         statsLenBeamCutouts[i][j] = statsTools.stats(clsloc)
    

# #         clsloc = [power_len_beam_cutouts[i][k][power] for k in arange(iMin, smallNum)]

# #         statsLenBeamCutouts[i][j] = statsTools.stats(clsloc)
    
        
# # prefac = lbin * (lbin+1) / 2 / numpy.pi



# # figure(1)
# # clf()
# # for j, power in enumerate(powers):
# #     subplot(nPowers, 1, j+1)
# #     for i, patch in enumerate(nsNames):

# #         errorbar(lbin + i * 10, prefac * statsLenBeamCutouts[i][j]['mean'], prefac * statsLenBeamCutouts[i][j]['stdev'] / sqrt(smallNum), label = patchNames[i], fmt = '.')

# #     plot(l_len, l_len**2 * theoryLensedPowers[j] / (2 * pi))

# #     xlim([0, 5000])
# #     ylabel('$l^2 ' + powersPretty[j] + r' / 2\pi$ (uK$^2$)')
# # legend()
# # show()
