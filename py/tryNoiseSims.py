

import liteMapPol
import sys
sys.path.append('/global/homes/e/engelen/aveTools/')
import aveTools
import cmblens.flipper.liteMap as liteMap
import numpy as np
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
    

import simTools
import matplotlib.pyplot as plt
# import cmblens.power


useQuickPower = True
doNoise = True
doCMB = False

doAll = True

if doCMB:
    print 'getting cmb sim'
    cmbSim = simTools.getActpolSim(iterationNum = 0, simType = 'cmb', doBeam = True)

    print cmbSim[0].data.shape





    cmbSimPower = aveTools.quickPower(cmbSim[0], \
                                      window = liteMapPol.initializeCosineWindow(cmbSim[0],200,0))
    cmbMap = liteMap.liteMapFromFits('/global/cscratch1/sd/engelen/s16maps/dataMaps/actpolDeep/dataCoadd_T_s13-d5-pa1-f150.fits')

    cmbPower = aveTools.quickPower(cmbMap, \
                               window = liteMapPol.initializeCosineWindow(cmbMap,200,0))
# patchList = [['deep1_s13_pa1_f150'] ,\
#              ['deep5_s13_pa1_f150'],\
#              ['deep6_s13_pa1_f150'],\
             
#              ['deep56_s14_pa1_f150'],\
#              ['deep56_s14_pa2_f150'],\
#              ['deep56_s15_pa1_f150'],\
#              ['deep56_s15_pa2_f150']]

# patchList = [['deep1_s13_pa1_f150'] ]
patchList = [['deep56_s15_pa3_f090'] , \
             ['deep56_s15_pa3_f150']]




patchList = [item for sublist in patchList for item in sublist]
nPatches = len(patchList)

nTqu = 3

if doNoise and doAll:


    if doAll:
        noiseSims = aveTools.onedl(nPatches)

    noiseMaps = aveTools.twodl(nPatches, nTqu)
    noiseSimPowers = aveTools.twodl(nPatches, nTqu)
    noisePowers = aveTools.twodl(nPatches, nTqu)
    
    for pi, patch in enumerate(patchList):
        print 'getting noise sim'
        if doAll:
            
            noiseSims[pi] = simTools.getActpolSim(iterationNum = 0, \
                                                  simType = 'noise', \
                                                  patch = patch)


        print 'patch', patch, 'shape is ', noiseSims[pi][0].data.shape


        crosslinkMap = liteMap.liteMapFromFits('/global/cscratch1/sd/engelen/s16maps/dataMaps/actpolDeep//crosslinkMap_T_' + patch + '.fits')

        noisePowerWindow = crosslinkMap.copy()
        noisePowerWindow.data *= (liteMapPol.initializeCosineWindow(crosslinkMap,200,0)).data

        for tqui in range(3):

            noiseMaps[pi][tqui] = liteMap.liteMapFromFits('/global/cscratch1/sd/engelen/s16maps/dataMaps/actpolDeep//noiseMap_%s_' % 'TQU'[tqui] + patch + '.fits')
        

            if useQuickPower:

                noiseSimPowers[pi][tqui] = aveTools.quickPower(noiseSims[pi][tqui], \
                                                               window = noisePowerWindow)
    
                
                noisePowers[pi][tqui] = aveTools.quickPower(noiseMaps[pi][tqui], \
                                                            window = noisePowerWindow)

            else:
                noiseSimPowers[pi][tqui] \
                    = cmblens.power.get_power(noiseSims[pi][tqui], \
                                              bin_file = 'binningTest',
                                              window = noisePowerWindow)
                
                noisePowers[pi][tqui]  \
                    = cmblens.power.get_poweraveTools.quickPower(noiseMaps[pi][tqui], \
                                                                 binFile = None,
                                                                 window = noisePowerWindow)

    

plt.clf()

if doCMB:
    plt.semilogy(cmbPower['lbin'], cmbPower['dlbin'], label = 'cmb power')
    plt.semilogy(cmbSimPower['lbin'], cmbSimPower['dlbin'], label = 'cmb sim power')
if doNoise:



    for pi, patch in enumerate(patchList):
        plt.figure(patch + ' ', figsize = (15,30))
        plt.clf()


        plt.subplot(4,2, 1)
        for tqui in range(3):
            plt.semilogy(noisePowers[pi][tqui]['lbin'], \
                         noisePowers[pi][tqui]['dlbin'], \
                         label = 'noise power %s' % 'TQU'[tqui],\
                         linestyle = 'dashed', \
                         color = 'grb'[tqui], linewidth = 4)
        
            plt.semilogy(noiseSimPowers[pi][tqui]['lbin'], \
                         noiseSimPowers[pi][tqui]['dlbin'], \
                         label = 'noise sim power %s' % 'TQU'[tqui], \
                         linestyle = 'solid', \
                         color = 'grb'[tqui], linewidth = 1.5)
        plt.legend(loc = 'best')
        plt.xlim([0,10000])
        plt.xlabel('$\ell$', fontsize = 20)
        plt.ylabel('$D_\ell$', fontsize = 20)

        plt.subplot(4,2, 2)
        for tqui in range(3):
            plt.plot(noisePowers[pi][tqui]['lbin'], \
                         noiseSimPowers[pi][tqui]['dlbin'] / noisePowers[pi][tqui]['dlbin'] , \
                         label = 'noise sim power / noise power %s' % 'TQU'[tqui],\
                         color = 'grb'[tqui], linewidth = 3)



        plt.xlim([0,10000])
        plt.legend()

        plt.axhline(1, linestyle = 'dotted', color = 'k')

        plt.xlim([0,10000])
        plt.xlabel('$\ell$', fontsize = 20)
        plt.ylabel('ratio', fontsize = 20)

        plt.gcf().text(0.5, 0.95, patch, fontsize=40, horizontalalignment = 'center')

        crosslinkMap = liteMap.liteMapFromFits('/global/cscratch1/sd/engelen/s16maps/dataMaps/actpolDeep//crosslinkMap_T_' + patch + '.fits')

        clims = [300, 300, 300]
        if True:
            for tqui in range(3):

                plt.subplot(4, 2, 3 + 2 * tqui)
                plt.imshow(crosslinkMap.data / np.mean(crosslinkMap.data) \
                           * noiseSims[pi][tqui].data , \
                           clim = clims[tqui] * np.array([-1,1]), cmap = 'RdBu')

                frame1 = plt.gca()
                frame1.axes.xaxis.set_ticklabels([])
                frame1.axes.yaxis.set_ticklabels([])


                plt.title('noise sim %s' % 'TQU'[tqui])

            for tqui in range(3):

                plt.subplot(4, 2, 3 + 2*tqui + 1 )
                plt.imshow(crosslinkMap.data / np.mean(crosslinkMap.data) \
                           * noiseMaps[pi][tqui].data, \
                           clim = clims[tqui] * np.array([-1,1]), cmap = 'RdBu')

                frame1 = plt.gca()
                frame1.axes.xaxis.set_ticklabels([])
                frame1.axes.yaxis.set_ticklabels([])
                
                
                plt.title('noise map %s' % 'TQU'[tqui])


        plt.savefig('../plot/trySimTools_%s.pdf' % patch)
stop

plt.figure(2)
plt.clf()

plt.semilogy(noisePower['lbin'], noisePower['dlbin'], label = 'noise power' )
             
plt.legend()
plt.show()

plt.figure(3)
plt.clf()
plt.imshow(cmbMap.data)

plt.show()


plt.figure(3)
plt.clf()
plt.imshow(noiseSim[0].data)
plt.show()
