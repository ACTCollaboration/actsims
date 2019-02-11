
#patch can be currently any of:
#                   [['deep1_s13_pa1_f150'] ,\
#                    ['deep5_s13_pa1_f150'],\
#                    ['deep6_s13_pa1_f150'],\
                   
#                    ['deep56_s14_pa1_f150'],\
#                    ['deep56_s14_pa2_f150'],\
#                    ['deep56_s15_pa1_f150'],\
#                    ['deep56_s15_pa2_f150'],\
#                    ['deep56_s15_pa3_f090'],\
#                    ['deep56_s15_pa3_f150'],\
#                    ['deep8_s15_pa1_f150'],\
#                    ['deep8_s15_pa2_f150'],\
#                    ['deep8_s15_pa3_f090'],\
#                    ['deep8_s15_pa3_f150']]


import orphics

import scipy
from actsims import simTools
from enlib import enmap
# import cmblens.power
import sys
sys.path.append("../../aveTools/")

import aveTools
import numpy as np

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import mapTools
import matplotlib.colors as colors
import pitas

import time

plt.switch_backend('Agg')

#Alternatively, you can pass the info separately, e.g: season = 's13', pa = 'pa1', freqGHz = 150, region = 'deep5'


                   # 'deep6_s13_pa1',\
                   #cut this one due to error in dict file


# patch = 'boss_s15_pa1_f150'
if True:
    psaList =     ['deep1_s13_pa1' ,\
                   'deep5_s13_pa1',\
                   'deep6_s13_pa1',\

                   'deep56_s14_pa1',\
                   'deep56_s14_pa2',\
                   
                   'deep56_s15_pa1',\
                   'deep56_s15_pa2',\
                   'deep56_s15_pa3',\

                   'deep8_s15_pa1',\
                   'deep8_s15_pa2',\
                   'deep8_s15_pa3',\

                   'boss_s15_pa1', \
                   'boss_s15_pa2', \
                   'boss_s15_pa3']
if False:
    # psaList =     ['deep1_s13_pa1']
        # psaList = ['deep56_s15_pa3']#  ,\
    psaList = ['boss_s15_pa3']#  ,\

    # psaList =     ['deep1_s13_pa1']

                   # 'deep5_s13_pa1',\
                   # 'deep6_s13_pa1']
if False:
    # psaList = ['deep56_s14_pa2',
    psaList = [ 'deep56_s15_pa2',\
               'deep56_s15_pa3',]

if False:
    # psaList = ['deep56_s14_pa2',
    psaList = ['deep56_s15_pa3',]
    # psaList = ['deep5_s13_pa1',]



import flipper.flipperDict as flipperDict 
p = flipperDict.flipperDict()
p.read_from_file('../inputParams/' + sys.argv[1])

mapSetArg = 'Mr3b'

overwrite = False

def my_smoothing(inmap):

    smoothed = scipy.ndimage.filters.gaussian_filter(inmap, 100)
    
    new = np.zeros(smoothed.shape)
    new[smoothed > .9] = 1.

    smoothedagain = scipy.ndimage.filters.gaussian_filter(new, 100)

    return smoothedagain

def kxKyFilter(inmap, kxcut = 90, kycut = 50, lmin = 0, lmax = 20000):
    
    from orphics import maps
    kmask = maps.mask_kspace(inmap.shape, inmap.wcs, kxcut, kycut, lmin, lmax)
    omap = maps.filter_map(inmap, kmask)
    return omap


import itertools

nMaps = len(psaList)

mapVersion = 'mr3c'


sDict = flipperDict.flipperDict()

sDict.read_from_file('../inputParams/signal.dict')
# sDict.read_from_file('../inputParams/templateInputs.dict')


beam150 = np.loadtxt(sDict['beamNames'][psaList[0] + '_' + 'f150'])

maxFreqs = 2
TQUs = 'TQU'
nTQUs = len(TQUs)

numIterations = 10
doFiltering = False
useSteveMasks = False


doAll = True
if doAll:
    simMaps = aveTools.onedl(len(psaList))
    simPowers = aveTools.fourdl(2, len(psaList), 3, numIterations)
    nullPowers = aveTools.threedl(2, len(psaList), 3)
    nullMaps = aveTools.threedl(2, len(psaList), 3)


    lmax      = 5000
    bin_edges = pitas.util.get_default_bin_edges(lmax)

    
    simCrossPowers = aveTools.fourdl(len(psaList) , maxFreqs * nTQUs , maxFreqs * nTQUs , numIterations )
    nullCrossPowers = aveTools.threedl(len(psaList) , maxFreqs * nTQUs , maxFreqs * nTQUs)



    freqs = [None] * len(psaList)
    for psai, psa in enumerate(psaList):
        for iii in range(numIterations):
        # for iii in [0]:

            print psa, iii
            # crossLinkMap = enmap.read_fits(p['crossLinkDict'][psa])

            # myPowerWindow  = my_smoothing(crossLinkMap)
            # myPowerWindow = crossLinkMap


            freqs[psai] = simTools.freqsInPsas(psa, p['freqsInArrays'])
            if useSteveMasks:
                steveDir = '/global/homes/s/schoi/projects/actpol/shared/for_alex/'
                myPowerWindow =  enmap.read_fits(steveDir + 'window/%s_f150_mr3_enki_ninkasi_win2_w0_cl0.00nK_pt1.00_nt0.0_T.fits' % psa)
            else:
                myPowerWindow = enmap.read_fits(p['crossLinkDict'][psa])
 # = enmap.read_fits(p['crossLinkDict'][psaPrime])

                if p['doWindowSmoothing']:
                    from actsims import noisePrep
                    if 'smoothingWidth' in p:
                        myPowerWindow = noisePrep.mySmoothing(myPowerWindow, N = p['smoothingWidth'])

                    else:
                        myPowerWindow = noisePrep.mySmoothing(myPowerWindow, N = p['smoothingWidth'])


            
            print 'warning assuming mask is same at all frequencies'
            maskName = 'mask_%s_%s_%s' %(psa, freqs[psai][0] ,
                                  ('stevewin' if useSteveMasks else 'smoothWin%d'%p['smoothingWidth']))
            print 'mask name is ', maskName
            pitas_lib = pitas.power.PITAS(maskName,
                myPowerWindow, myPowerWindow, bin_edges, lmax, None, overwrite)




            for freqi, freq in enumerate(freqs[psai]):
                for tqui, tqu in enumerate(TQUs):
                    nullMaps[freqi][psai][tqui] \
                        = enmap.read_map(p['dataMapDir'] + 'nullMaps' \
                                         + 'IQU'[tqui] + '_' + psa + '_' + freqs[psai][freqi] + '.fits')

                    if doFiltering:
                        nullMaps[freqi][psai][tqui] = kxKyFilter(np.nan_to_num(nullMaps[freqi][psai][tqui]))

            print 'getting noise sim for  ', psa, 'noiseDictFile is ', sys.argv[1]
            start = time.clock()
            simMaps[psai] = simTools.getActpolSim(iterationNum = iii + 1, \
                                                  simType = 'noise', \
                                                  psa = psa)  #, coaddDictFile = 'Coadd_boss_test.dict')


                    # simToolsDictFile = 'simToolsPrime.dict')
                # simMaps[mt] = enmap.ndmap(np.zeros(crossLinkMap.shape), crossLinkMap.wcs)
            print 'done, took' , time.clock() - start , 'seconds'

            if doFiltering:
                simMaps[psai] = kxKyFilter(simMaps[psai])

            print 'getting all the cross spectra for ' , psa ,'and iteration ', iii
            start = time.clock()
            for index, (freq, tqu) in enumerate(itertools.product(freqs[psai], TQUs)):
                for indexPrime, (freqPrime, tquPrime) in enumerate(itertools.product(freqs[psai], TQUs)):
                    # if index != indexPrime :
                    #     continue
                    
                    freqi = freqs[psai].index(freq)
                    freqiPrime = freqs[psai].index(freqPrime)
                    tqui = TQUs.index(tqu)
                    tquiPrime = TQUs.index(tquPrime)

                    print 'doing ' , freqi, freqiPrime, tqui, tquiPrime


                    simCrossPowers[psai][index][indexPrime][iii] = {}

                    simCrossPowers[psai][index][indexPrime][iii]['lbin'], \
                        simCrossPowers[psai][index][indexPrime][iii]['clbin'] \
                        = pitas_lib.get_power(
                            emap1 = np.nan_to_num(np.squeeze(simMaps[psai][freqi, tqui, :, :]) \
                                                  * myPowerWindow ),
                            emap2 = np.nan_to_num(np.squeeze(simMaps[psai][freqiPrime, tquiPrime, :, :]) \
                                                  * myPowerWindow))

                    nullCrossPowers[psai][index][indexPrime] = {}

                    nullCrossPowers[psai][index][indexPrime]['lbin'], \
                        nullCrossPowers[psai][index][indexPrime]['clbin']  \
                        = pitas_lib.get_power(
                            emap1 = np.nan_to_num(np.squeeze(nullMaps[freqi][psai][tqui]) \
                                                  * myPowerWindow ),
                            emap2 = np.nan_to_num(np.squeeze(nullMaps[freqiPrime][psai][tquiPrime])  \
                            * myPowerWindow))


                    if freq == freqPrime and tqu == tquPrime:


                        simPowers[freqi][psai][tqui][iii] = simCrossPowers[psai][index][indexPrime][iii]
                        nullPowers[freqi][psai][tqui] = nullCrossPowers[psai][index][indexPrime]
            print 'done, took' , time.clock() - start, 'seconds'
    #now delete the maps, they are taking up too much space
    
    # nullMaps[psai] = None
    # simMaps[psai] = None
                # nullMapOld = enmap.read_map('/global/cscratch1/sd/engelen/s16maps/dataMaps/actpolDeep/nullMapsI_deep56_s14_pa2_f150.fits')
                # nullMapOldPower = aveTools.quickPowerEnlib(nullMapOld , window = myPowerWindow)

                # simoneSim =  enmap.read_map('/global/cscratch1/sd/engelen/noise_hack/noise_sims_test/s14_mr2_deep56_pa2_f150_night_noise_test_20180417.fits')
                # simonePower = aveTools.quickPowerEnlib(simoneSim , window = myPowerWindow)

        # nullMapCutout  = enmap.read_fits('/global/cscratch1/sd/engelen/s16maps/dataMaps/actpolDeep/nullMapsI_deep56_s14_pa2_f150.fits', box = box)

if True:
    for psai, psa in enumerate(psaList):
        plt.figure(psa + 'all crosses', figsize = (8,8))
        j = 0
        for index, (freq, tqu) in enumerate(itertools.product(freqs[psai], TQUs)):
            for indexPrime, (freqPrime, tquPrime) in enumerate(itertools.product(freqs[psai], TQUs)):

                j += 1
                if index > indexPrime:
                    continue
                lbin = nullPowers[freqi][psai][tqui]['lbin']
                fullStats = aveTools.stats(simCrossPowers[psai][index][indexPrime], label = 'clbin')

                nPanels = len(freqs[psai]) * len(TQUs)

                plt.subplot(nPanels,nPanels, j)

                plt.semilogy(lbin, fullStats['mean'], \
                             label = 'noise sim' ,\
                             linestyle = 'solid', \
                             linewidth = 1.5)

                plt.semilogy(nullCrossPowers[psai][index][indexPrime]['lbin'], \
                             nullCrossPowers[psai][index][indexPrime]['clbin'], \
                             label = 'data null',\
                             linestyle = 'dashed', \
                             linewidth = 4)
                plt.title( '%s%s X %s%s' % (tqu, freq,  tquPrime, freqPrime ))
                
                if j == 1:
                    plt.legend()
        plt.tight_layout()
        plt.savefig('../plot/%s_allcrosses.pdf' % psa)
if True:



    for psai, psa in enumerate(psaList):
        nPanels = len(freqs[psai]) * len(TQUs)
        fullStats = aveTools.twodl(nPanels, nPanels)        

        for index, (freq, tqu) in enumerate(itertools.product(freqs[psai], TQUs)):
            for indexPrime, (freqPrime, tquPrime) in enumerate(itertools.product(freqs[psai], TQUs)):
                fullStats[index][indexPrime] = aveTools.stats(simCrossPowers[psai][index][indexPrime], label = 'clbin')

        plt.figure(psa + 'correlations', figsize = (8,8))
        j = 0
        for index, (freq, tqu) in enumerate(itertools.product(freqs[psai], TQUs)):
            for indexPrime, (freqPrime, tquPrime) in enumerate(itertools.product(freqs[psai], TQUs)):

                j += 1
                if index >= indexPrime:
                    continue
                lbin = nullPowers[freqi][psai][tqui]['lbin']

                plt.subplot(nPanels,nPanels, j)

                plt.plot(lbin, fullStats[index][indexPrime]['mean'] \
                         / np.sqrt(fullStats[index][index]['mean'] * fullStats[indexPrime][indexPrime]['mean']), \
                         label = 'noise sim',\
                         linestyle = 'solid', \
                         linewidth = 1.5)

                plt.plot(lbin,  \
                             nullCrossPowers[psai][index][indexPrime]['clbin'] \
                             / np.sqrt(nullCrossPowers[psai][index][index]['clbin'] * nullCrossPowers[psai][indexPrime][indexPrime]['clbin']), 
                             label = 'data null',\
                             linestyle = 'dashed', \
                             linewidth = 4)

                plt.title( '%s%s X %s%s' % (tqu, freq,  tquPrime, freqPrime ))
                plt.axhline(0, linestyle = 'dotted', linewidth = .3, color = 'k')
                if j == 2:
                    plt.legend()

                plt.ylim([-.6,.6])
        plt.tight_layout()
        plt.savefig('../plot/%s_correlations.pdf' % psa)
        

                
                

if True:
    for psai, psa in enumerate(psaList):

        for freqi, freq in enumerate(freqs[psai]):
            print 'psa, freq', psa, freq
            plt.figure(psa, figsize = (15,8))
            plt.clf()

            plt.subplot(1,2, 1)

            for tqui in range(3):
                lbin = nullPowers[freqi][psai][tqui]['lbin']
                theseStats = aveTools.stats(simPowers[freqi][psai][tqui], label = 'clbin')

                plt.semilogy(lbin, \
                             theseStats['mean'], \
                             label = 'noise sim   %s %s' % ('TQU'[tqui], freqs[psai][freqi]),\
                             linestyle = 'solid', \
                             color = 'grb'[tqui], linewidth = 1.5)
                plt.semilogy(nullPowers[freqi][psai][tqui]['lbin'], \
                             nullPowers[freqi][psai][tqui]['clbin'], \
                             label = 'data null  %s %s' % ('TQU'[tqui], freqs[psai][freqi]),\
                             linestyle = 'dashed', \
                             color = 'grb'[tqui], linewidth = 4)

            # plt.semilogy(nullMapOldPower['lbin'], nullMapOldPower['clbin'], label = 'nullMapOld')
            # plt.semilogy(simonePower['lbin'], simonePower['clbin'], label = 'simone sim')


            plt.legend(loc = 'best')
            plt.xlim([0,10000])
            plt.xlabel('$\ell$', fontsize = 20)
            plt.ylabel('$C_\ell$', fontsize = 20)

            plt.subplot(1, 2, 2)

            for tqui in range(3):
                theseStats = aveTools.stats(simPowers[freqi][psai][tqui], label = 'clbin')

                plt.errorbar(lbin + tqui / 3. * (lbin[1] - lbin[0]), \
                             theseStats['mean'] / nullPowers[freqi][psai][tqui]['clbin'], \
                             theseStats['stdev'] / np.sqrt(numIterations) / nullPowers[freqi][psai][tqui]['clbin'], \

                             label = '(noise sim) / (data null)   %s %s' % ('TQU'[tqui], freqs[psai][freqi]),\
                             # linestyle = 'solid', \
                             color = 'grb'[tqui], linewidth = 1.5, capsize = 0 , linestyle = 'none', marker = 'o', markersize = 2.5, markeredgewidth = 0)


            plt.xlim([0,10000])
            plt.legend()

            plt.axhline(1, linestyle = 'dotted', color = 'k')

            plt.xlim([0,10000])
            plt.xlabel('$\ell$', fontsize = 20)
            plt.ylabel('ratio', fontsize = 20)

            plt.gcf().text(0.5, 0.95, psa + ' ' + freqs[psai][freqi] + ' ' + mapVersion, fontsize=40, horizontalalignment = 'center')

            plt.ylim([.5, 1.5])
            plt.savefig('../plot/sim_vs_nullmaps_%s_%s_%s.pdf' %(psa, freq, mapVersion))

            # plt.show()
stop            
lmax = 10000


plt.figure('steve vs alex')
plt.clf()

tqui = 0
for psai, psa in enumerate(psaList):

    for freqi, freq in enumerate(freqs[psai]):

        steveSimPower = np.loadtxt('../data/%s_%s_sim_noise_ps_180728.txt' % (psa, freq))

        plt.subplot(2, 2, 1)

            
        plt.semilogy(steveSimPower[:, 0], steveSimPower[:, 1] \
                     / (steveSimPower[:, 0] * (steveSimPower[:, 0] + 1) / (2. * np.pi)), label = 'steve sim power')

        theseStats = aveTools.stats(simPowers[freqi][psai][tqui], label = 'clbin')

        plt.semilogy(lbin, \
                     theseStats['mean'], \
                     label = 'alex  sim   %s %s' % ('TQU'[tqui], freqs[psai][freqi]),\
                     color = 'grb'[tqui], linewidth = 4, linestyle = 'dashed')

        plt.legend()
        plt.ylabel('$C_\ell$')
        steveNullPower = np.loadtxt('../data/%s_%s_data_noise_ps_180728.txt' % (psa, freq))



        plt.subplot(2, 2, 3)

        alexSteveSimRatio = theseStats['mean'] /  \
                            scipy.interpolate.interp1d(steveSimPower[:, 0],  steveSimPower[:, 1] \
                                                       / (steveSimPower[:, 0] * (steveSimPower[:, 0] + 1) / (2. * np.pi)),
                                                       bounds_error = False )(lbin)
                 
            
        plt.plot(lbin, alexSteveSimRatio ,
                 label = 'alex / steve, sims')
        plt.plot(beam150[:, 0], beam150[:, 1]**2, label = 'beam squared')
        plt.xlim([0,8e3])


        plt.legend()
        plt.ylabel('relative')
        

        plt.subplot(2, 2, 2)

            
        plt.semilogy(steveNullPower[:, 0], steveNullPower[:, 1] \
                     / (steveNullPower[:, 0] * (steveNullPower[:, 0] + 1) / (2. * np.pi)), label = 'steve data-noise power')

        plt.semilogy(nullPowers[freqi][psai][tqui]['lbin'], \
                     nullPowers[freqi][psai][tqui]['clbin'], \
                     label = 'alex  %s %s' % ('TQU'[tqui], freqs[psai][freqi]),\
                     linestyle = 'dashed', \
                     color = 'grb'[tqui],  linewidth = 4)

        plt.legend()
        plt.ylabel('$C_\ell$')

        plt.subplot(2, 2, 4)

        alexSteveNullRatio = nullPowers[freqi][psai][tqui]['clbin'] / \
                             scipy.interpolate.interp1d(steveNullPower[:, 0],  steveNullPower[:, 1] \
                                                        / (steveNullPower[:, 0] * (steveNullPower[:, 0] + 1) / (2. * np.pi)) ,
                                                        bounds_error = False )(lbin)

        plt.plot(lbin,alexSteveNullRatio,
                 label = 'alex / steve, data-noise')

        plt.xlim([0,8e3])

        plt.plot(beam150[:, 0], beam150[:, 1]**2, label = 'beam squared')

        plt.legend()

stop
plt.show()

plt.figure('test')
plt.clf()
plt.plot(lbin, alexSteveSimRatio / alexSteveNullRatio)
plt.ylim([.8, 1.2])
plt.show()

stop


plotTypes = ['nullMap', 'simMap']



for pti, plotType in enumerate(plotTypes):

    for psai, psa in enumerate(psaList):
        # if pti == 0:
        #     climVals = np.zeros(3 * simMaps[psai].shape[0] * 
        #                         3 * simMaps[psai].shape[0] + 1) #+1 due to the fact that j starts at 1

        j = 1

        plt.figure(psa + ' 2d sim', figsize = (15,15))
        # plt.figure(psa + ' 2d data', figsize = (15,15))

        plt.clf()

        for index, (freq, tqu) in enumerate(itertools.product(freqs[psai], TQUs)):
            for indexPrime, (freqPrime, tquPrime) in enumerate(itertools.product(freqs[psai], TQUs)):
                freqi = freqs[psai].index(freq)
                freqiPrime = freqs[psai].index(freqPrime)
                tqui = TQUs.index(tqu)
                tquiPrime = TQUs.index(tquPrime)

                if indexPrime > index:
                    j+= 1
                    continue
                

                plt.subplot(2, 1, j + 1)# 3 * simMaps[psai].shape[0],
                            # 3 * simMaps[psai].shape[0],
                            # j
                # )

                # if indexPrime > index:
                #     j += 1

                #     continue

                # print 'doing ' , patch , 'with ', patchPrime, '; doing ', iqu , 'with ', iquPrime


                # crossLinkMap = enmap.read_fits(p['crossLinkDict'][psa])

                # fc = orphics.maps.FourierCalc(simMaps[psai].shape, simMaps[psai].wcs)
                # power2d, kmap1, kmap2 =  fc.power2d(np.nan_to_num(simMaps[psai][freqi, tqui, :, :]),
                #                                     np.nan_to_num(simMaps[psai][freqiPrime, tquiPrime, :, :]))
                if plotType == 'nullMap':
                    power2d = nullCrossPowers[psai][index][indexPrime]['power2d']
                else:
                    power2d = simCrossPowers[psai][index][indexPrime]['power2d']

                lmap = power2d.modlmap()

                lx = lmap[1,0]
                ly = lmap[0,1]

                # toPlot = power2d[0 : lmax / lx, -lmax / ly:lmax / ly]
                toPlot = power2d

                #only set clim for first pass.
                if pti == 0:
                    climVals[j] =  np.max([toPlot.max(), -toPlot.max()])
                plt.imshow(np.flipud((((np.fft.fftshift(toPlot))))) , extent = [-lmap[:,0].max(), lmap[:,0].max(),
                                                                                -lmap[0,:].max(), lmap[0,:].max()],
                           interpolation = 'none', norm = colors.SymLogNorm(linthresh = climVals[j] / 1e5),
                           cmap = 'RdBu',
                           clim = [-climVals[j], +climVals[j]])
                # plt.imshow(np.flipud(((np.log10(np.fft.fftshift(power2d))))) )
                plt.colorbar()
                plt.title('%s %s  X  %s %s' % (freq, tqu, freqPrime, tquPrime), fontsize = 10)

                j += 1

        plt.gcf().text(0.5, 0.95, psa + ' ' + plotType, fontsize=40, horizontalalignment = 'center')

        plt.savefig('../plot/twoDSpectra_sims_' + psa + '_' + plotType + '.png')

for psai, psa in enumerate(psaList):
    j = 1

    plt.figure(psa + ' correlations', figsize = (15,15))
    # plt.figure(psa + ' 2d data', figsize = (15,15))

    plt.clf()

    for index, (freq, tqu) in enumerate(itertools.product(freqs[psai], TQUs)):
        for indexPrime, (freqPrime, tquPrime) in enumerate(itertools.product(freqs[psai], TQUs)):
            freqi = freqs[psai].index(freq)
            freqiPrime = freqs[psai].index(freqPrime)
            tqui = TQUs.index(tqu)
            tquiPrime = TQUs.index(tquPrime)
            
            if indexPrime >= index:
                j += 1
                continue

            plt.subplot(3 * simMaps[psai].shape[0],
                        3 * simMaps[psai].shape[0],
                        j)

            for  plotType in plotTypes:
                if plotType == 'nullMap':
                    powerDict = nullCrossPowers
                else:
                    powerDict = simCrossPowers

                plt.plot(powerDict[psai][index][indexPrime]['lbin'],
                              powerDict[psai][index][indexPrime]['clbin'] \
                              / np.sqrt(powerDict[psai][index][index]['clbin'] * powerDict[psai][indexPrime][indexPrime]['clbin']),
                             label = plotType)
                plt.axhline(linestyle =  'dotted')
                plt.title('%s %s  X  %s %s' % (freq, tqu, freqPrime, tquPrime), fontsize = 10)
                plt.ylim([-.4, .4])
                # plt.yscale('symlog')
            j += 1
    plt.legend()
    
    plt.gcf().text(0.5, 0.95, psa, fontsize=40, horizontalalignment = 'center')

    plt.gcf().subplots_adjust(left=None, bottom=None, right=None, top=None,
                                    wspace=.3, hspace=.3)
    plt.savefig('../plot/correlations_' + psa + '.pdf')


stop

colors = ['r', 'b', 'g']
plt.figure('hack sim powers')
plt.clf()
plt.subplot(2,1,1)


for mt, mapType in enumerate(mapTypes):

    if mapType == 'mat':
        continue

    plt.semilogy(simPowersMCM[mt]['lbin'], simPowersMCM[mt]['clbin'], label = mapType, color = colors[mt])

plt.semilogy(nullPowerMCM['lbin'], nullPowerMCM['clbin'], label = 'data null')

plt.legend()
plt.ylabel('$C_\ell$')

plt.subplot(2,1,2)

for mt, mapType in enumerate(mapTypes):
    if mapType == 'mat':
        continue
    plt.plot(simPowersMCM[mt]['lbin'], simPowersMCM[mt]['clbin'] / nullPowerMCM['clbin'] , label = mapType + ' / (data null)', color = colors[mt])



    
plt.legend()
plt.axhline(1., color = '.5', linestyle = 'dashed', linewidth = 2)
# plt.show()
plt.savefig('../plot/simone_mat_alex.png')

            # clims = [300, 300, 300]

            # for tqui in range(3):

            #     plt.subplot(4, 2, 3 + 2 * tqui)
            #     plt.imshow( np.squeeze(simMaps[psai][freqi, tqui, :, :]) , \
            #                clim = clims[tqui] * np.array([-1,1]), cmap = 'RdBu')

            #     frame1 = plt.gca()
            #     frame1.axes.xaxis.set_ticklabels([])
            #     frame1.axes.yaxis.set_ticklabels([])

            #     plt.title('noise sim %s' % 'TQU'[tqui])

            # for tqui in range(3):

            #     plt.subplot(4, 2, 3 + 2*tqui + 1 )
            #     plt.imshow( np.squeeze(nullMaps[freqi][psai][ tqui]) , \
            #                clim = clims[tqui] * np.array([-1,1]), cmap = 'RdBu')


            #     frame1 = plt.gca()
            #     frame1.axes.xaxis.set_ticklabels([])
            #     frame1.axes.yaxis.set_ticklabels([])


            #     plt.title('data null %s' % 'TQU'[tqui])
