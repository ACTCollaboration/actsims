

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

doCmb = True

doAll = True

patchList = [['deep1_s13_pa1_f150'] ,\
             ['deep5_s13_pa1_f150'],\
             ['deep6_s13_pa1_f150'],\
             
             ['deep56_s14_pa1_f150'],\
             ['deep56_s14_pa2_f150'],\
             ['deep56_s15_pa1_f150'],\
             ['deep56_s15_pa2_f150']]

# # patchList = [['deep1_s13_pa1_f150'] ]
# patchList = [['deep56_s15_pa3_f090'] , \
#              ['deep56_s15_pa3_f150']]
# patchList = [['deep1_s13_pa1_f150']]

patchList = [item for sublist in patchList for item in sublist]
nPatches = len(patchList)

nTqu = 3
if doCmb and doAll:


    if doAll:
        cmbSims = aveTools.onedl(nPatches)

    cmbSimPowers = aveTools.onedl(nPatches)

    
    for pi, patch in enumerate(patchList):
        print 'getting cmb sim'
        if doAll:
            
            cmbSims[pi] = simTools.getActpolSim(iterationNum = 0, \
                                                simType = 'cmb', \
                                                patch = patch, \
                                                doBeam = True)


        print 'patch', patch, 'shape is ', cmbSims[pi][0].data.shape


        # crosslinkMap = liteMap.liteMapFromFits('/global/cscratch1/sd/engelen/s16maps/dataMaps/actpolDeep//crosslinkMap_T_' + patch + '.fits')

        # noisePowerWindow = crosslinkMap.copy()
        # noisePowerWindow.data *= (liteMapPol.initializeCosineWindow(crosslinkMap,200,0)).data
        cmbPowerWindow = cmbSims[pi][0].copy()
        cmbPowerWindow.data = (liteMapPol.initializeCosineWindow(cmbSims[pi][0],200,0)).data

        if useQuickPower:
            cmbSims[pi][1].data *= -1.
            cmbSimPowers[pi] = aveTools.allpowers(cmbSims[pi][0], \
                                                  cmbSims[pi][1], \
                                                  cmbSims[pi][2], \
                                                  window = cmbPowerWindow, \
                                                  binFile = 'binningTest')
    
        else:
            raise ValueError('not yet implemented')




if doCmb:

    specs = ['cl_TT', 'cl_EE', 'cl_BB']
    plt.figure('cmb powers', figsize = (8,8))
    plt.clf()
    for pi, patch in enumerate(patchList):

        for si, spec in enumerate(specs):
            # plt.subplot(3,1,si + 1)
            plt.semilogy(cmbSimPowers[pi]['lbin'], \
                         cmbSimPowers[pi][spec] \
                         * cmbSimPowers[pi]['lbin'] \
                         * (cmbSimPowers[pi]['lbin']  + 1) / (2. * np.pi), \
                         label = patch + ' ' + spec if si == 0 else None,\
                         # linestyle = 'dashed', \
                         color = 'grb'[si], linewidth = .5)
        
    plt.legend(loc = 'best')
    plt.xlim([0,10000])
    plt.ylim([1e-7,1e4])

    plt.xlabel('$\ell$', fontsize = 20)
    plt.ylabel('$D_\ell$', fontsize = 20)


    plt.savefig('../plot/trySignalSims_%s.pdf' % patch)
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
