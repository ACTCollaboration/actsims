import numpy as np
import scipy
import matplotlib.pyplot as plt
from soapack import interfaces
from pixell import enmap
from pixell import enplot
import os,sys
from orphics import mpi,io
from pitas import power as ppower
from enlib import bench

"""
This script calculates all ACT and Planck cross-spectra.
"""

import argparse
# Parse command line
parser = argparse.ArgumentParser(description='Do a thing.')
parser.add_argument("region", type=str,help='boss or deep56?')
parser.add_argument("--fft", action='store_true',help='A flag.')
args = parser.parse_args()

fftstr = "_fft" if args.fft else ""

mask_version = interfaces.dconfig['act_mr3']['default_mask_version']
mask = interfaces.get_act_mr3_crosslinked_mask(args.region)
identifier = args.region+"_"+mask_version

comm = mpi.MPI.COMM_WORLD
rank = comm.Get_rank()
numcores = comm.Get_size()

#bin_edges = np.arange(80,8000,100)
lefts,rights = np.loadtxt("data/BIN_ACTPOL_50_4_SC_low_ell_startAt2.txt",unpack=True,usecols=[0,1])
bin_edges = np.append(lefts,rights[-1:])
print(bin_edges)


# Calculate and store mcm
if not(args.fft):
    with io.nostdout():
        if rank==0:
            pcalc = ppower.PITAS(identifier,mask,mask,bin_edges)
        comm.Barrier()
        pcalc = ppower.PITAS(identifier,mask,mask,bin_edges)



# We seek to plot some ACTxACT, PlanckxPlanck, and ACTxPlanck cross spectra. Here we get the coadds from the two Planck half-missions.
def get_planck_coadd(freq, dmp):
    psplit0 = dmp.get_split(freq, splitnum=0, ncomp=1, srcfree=True)
    psplit1 = dmp.get_split(freq, splitnum=1, ncomp=1, srcfree=True)
    psplit0_i = dmp.get_split_ivar(freq, splitnum=0, ncomp=1, srcfree=True)
    psplit1_i = dmp.get_split_ivar(freq, splitnum=1, ncomp=1, srcfree=True)
    weighted = (psplit0_i * psplit0 + psplit1 * psplit1_i) / (psplit0_i + psplit1_i)
    weighted[~np.isfinite(weighted)] = 0.0
    return weighted


# Here is where we do a crappy job of estimating power spectra.


def bin(data,modlmap,bin_edges):
    digitized = np.digitize(np.ndarray.flatten(modlmap), bin_edges,right=True)
    return np.bincount(digitized,(data).reshape(-1))[1:-1]/np.bincount(digitized)[1:-1]


def compute_ps(map1, map2, beamf1, beamf2):
    """Compute the FFTs, multiply, bin
    """
    if args.fft:
        kmap1 = enmap.fft(map1*mask, normalize="phys")
        kmap2 = enmap.fft(map2*mask, normalize="phys")
        power = (kmap1*np.conj(kmap2)).real
        bin_edges = np.arange(20,8000,40)
        centers = (bin_edges[1:] + bin_edges[:-1])/2.
        w2 = np.mean(mask**2.)
        modlmap = enmap.modlmap(map1.shape,map1.wcs)
        binned_power = bin(power/w2/beamf1(modlmap)/beamf2(modlmap),modlmap,bin_edges)
        return centers, binned_power
    else:
        ells,cls = pcalc.get_power_scalarXscalar(map1*mask, map2*mask,ret_dl=False)
        return ells,cls/beamf1(ells)/beamf2(ells)


# # ACTxPlanck

# coadded ACT x coadded Planck


patch = args.region
dma = interfaces.ACTmr3(region=mask)
dmp = interfaces.PlanckHybrid(region=mask)
pfreqs = ['030','044','070','100','143','217','353','545']
nfreqs = len(pfreqs)

# we loop over all pairs of Planck x ACT

# combs = []
# for planckfreq in ['030','044','070','100','143','217','353','545']: # no '857'
#     for actseason in ['s14','s15']:
#         for array in ['pa1_f150', 'pa2_f150', 'pa3_f090', 'pa3_f150']:
#             try:
#                 actbeam = lambda x: dma.get_beam(x, actseason, 
#                                        patch, array)
#                 actbeam(100)
#                 combs.append((planckfreq,actseason,array))
#             except:
#                 pass

# Njobs = len(combs)
# num_each,each_tasks = mpi.mpi_distribute(Njobs,numcores,allow_empty=True)
# if rank==0: print ("At most ", max(num_each) , " tasks...")
# my_tasks = each_tasks[rank]

# for task in my_tasks:
#     planckfreq,actseason,array = combs[task]
#     print(planckfreq,actseason,array)
#     planckbeam = lambda x: dmp.get_beam(x, planckfreq)
#     planckmap = get_planck_coadd(planckfreq, dmp)[0,:,:]
#     actbeam = lambda x: dma.get_beam(x, actseason, 
#                            patch, array)
#     actmap = dma.get_coadd(actseason, 
#                            patch, array, ncomp=1, 
#                            srcfree=True)[0,:,:] # just want T
#     lb, Cb = compute_ps(planckmap, actmap, planckbeam, actbeam)
#     io.save_cols("/scratch/r/rbond/msyriac/data/depot/actsims/spectra/spec%s_%s_%s_%s_%s.txt" % (fftstr,patch,planckfreq,actseason,array),(lb,Cb))
#     if rank==0: print ("Rank 0 done with task ", task+1, " / " , len(my_tasks))

# comm.Barrier()
# # # Planck x Planck (different freqs)
# # I use coadded planck x coadded planck

# combs = []
# for i in range(nfreqs):
#     for j in range(i+1,nfreqs):
#         combs.append((i,j))

# Njobs = len(combs)
# num_each,each_tasks = mpi.mpi_distribute(Njobs,numcores,allow_empty=True)
# if rank==0: print ("At most ", max(num_each) , " tasks...")
# my_tasks = each_tasks[rank]

# for task in my_tasks:
#     i,j = combs[task]
#     planckfreq0 = pfreqs[i]
#     planckfreq1 = pfreqs[j]

#     planckbeam0 = lambda x: dmp.get_beam(x, planckfreq0)
#     planckmap0 = get_planck_coadd(planckfreq0, dmp)[0,:,:]
#     planckbeam1 = lambda x: dmp.get_beam(x, planckfreq1)
#     planckmap1 = get_planck_coadd(planckfreq1, dmp)[0,:,:]
#     lb, Cb = compute_ps(planckmap0, planckmap1, planckbeam0, planckbeam1)
#     io.save_cols("/scratch/r/rbond/msyriac/data/depot/actsims/spectra/spec%s_%s_%s_%s.txt" % (fftstr,patch,planckfreq0,planckfreq1),(lb,Cb))
#     if rank==0: print ("Rank 0 done with task ", task+1, " / " , len(my_tasks))


# comm.Barrier()
# # # Planck x Planck (same freq)
# # These are spectra for which $f_0 = f_1$, so I use half missions.


# Njobs = len(pfreqs)
# num_each,each_tasks = mpi.mpi_distribute(Njobs,numcores,allow_empty=True)
# if rank==0: print ("At most ", max(num_each) , " tasks...")
# my_tasks = each_tasks[rank]

# for task in my_tasks:
#     planckfreq = pfreqs[task]
#     # we loop over all pairs of Planck x Planck
#     planckbeam = lambda x: dmp.get_beam(x, planckfreq)
#     planckmap0 = dmp.get_split(planckfreq, splitnum=0, ncomp=1, srcfree=True)
#     planckmap1 = dmp.get_split(planckfreq, splitnum=1, ncomp=1, srcfree=True)

#     lb, Cb = compute_ps(planckmap0, planckmap1, planckbeam, planckbeam)

#     io.save_cols("/scratch/r/rbond/msyriac/data/depot/actsims/spectra/spec%s_%s_%s_%s.txt" % (fftstr,patch,planckfreq,planckfreq),(lb,Cb))
#     if rank==0: print ("Rank 0 done with task ", task+1, " / " , len(my_tasks))


# comm.Barrier()
# # ACT x ACT
# Different seasons/arrays - can use just the coadds


combs = []
    
# we loop over all pairs of ACT x ACT
for actseason0 in ['s14','s15']: # s13 doesn't have these patches
    for array0 in ['pa1_f150', 'pa2_f150', 'pa3_f090', 'pa3_f150']:
        for actseason1 in ['s14','s15']: # s13 doesn't have these patches
            for array1 in ['pa1_f150', 'pa2_f150', 'pa3_f090', 'pa3_f150']:
                try:
                    if (actseason0==actseason1) and (array0==array1): raise
                    if array0=='pa3_f150' and array1=='pa3_f090': raise
                    if array1=='pa3_f150' and array0=='pa3_f090': raise
                    actbeam0 = lambda x: dma.get_beam(x, actseason0, patch, array0)
                    actbeam1 = lambda x: dma.get_beam(x, actseason1, patch, array1)
                    actbeam0(100)
                    actbeam1(100)
                    print(actseason0,array0,actseason1,array1)
                    combs.append((actseason0,array0,actseason1,array1))
                except:
                    pass

Njobs = len(combs)
num_each,each_tasks = mpi.mpi_distribute(Njobs,numcores,allow_empty=True)
if rank==0: print ("At most ", max(num_each) , " tasks...")
my_tasks = each_tasks[rank]

for task in my_tasks:
    actseason0,array0,actseason1,array1 = combs[task]
    actbeam0 = lambda x: dma.get_beam(x, actseason0, patch, array0)
    actbeam1 = lambda x: dma.get_beam(x, actseason1, patch, array1)

    actmap0 = dma.get_coadd(actseason0, patch, array0, 
                            ncomp=1, srcfree=True)[0,:,:] # just want T
    actmap1 = dma.get_coadd(actseason1, patch, array1, 
                            ncomp=1, srcfree=True)[0,:,:] # just want T
    lb, Cb = compute_ps(actmap0, actmap1, actbeam0, actbeam1)

    io.save_cols("/scratch/r/rbond/msyriac/data/depot/actsims/spectra/spec%s_%s_%s_%s_%s_%s.txt" % (fftstr,patch, actseason0, array0, actseason1, array1),(lb,Cb))
    if rank==0: print ("Rank 0 done with task ", task+1, " / " , len(my_tasks))


comm.Barrier()

# combs = []
# nsplits = 4

# # we loop over all pairs of ACT x ACT
# for actseason0 in ['s14','s15']: # s13 doesn't have these patches
#     for array0 in ['pa1_f150', 'pa2_f150', 'pa3_f090', 'pa3_f150']:
#         try:
#             actbeam = lambda x: dma.get_beam(x, actseason0, patch, array0)
#             actbeam(100)
#             combs.append((actseason0,array0,actseason0,array0))
#         except:
#             pass
# combs.append(("s15","pa3_f090","s15","pa3_f150"))

        
# Njobs = len(combs)
# num_each,each_tasks = mpi.mpi_distribute(Njobs,numcores,allow_empty=True)
# if rank==0: print ("At most ", max(num_each) , " tasks...")
# my_tasks = each_tasks[rank]

# for task in my_tasks:
#     actseason0,array0,actseason1,array1 = combs[task]
            
#     actbeam0 = lambda x: dma.get_beam(x, actseason0, patch, array0)
#     actbeam1 = lambda x: dma.get_beam(x, actseason1, patch, array1)
#     actmaps0 = dma.get_splits(actseason0, patch, array0, 
#                                 ncomp=1, srcfree=True)
#     actmaps1 = dma.get_splits(actseason1, patch, array1, 
#                                 ncomp=1, srcfree=True)

#     spec = 0.
#     count = 0
#     for i in range(nsplits):
#         for j in range(i+1,nsplits):
#             actmap0 = actmaps0[0, i, 0, :, :]
#             actmap1 = actmaps1[0, j, 0, :, :]
#             lb, Cb = compute_ps(actmap0, actmap1, actbeam0, actbeam1)
#             spec = spec + np.asarray(Cb)
#             count = count + 1

#     Cb = spec / count
#     io.save_cols("/scratch/r/rbond/msyriac/data/depot/actsims/spectra/spec%s_%s_%s_%s_%s_%s.txt" % (fftstr,patch, actseason0, array0, actseason1, array1),(lb,Cb))
#     if rank==0: print ("Rank 0 done with task ", task+1, " / " , len(my_tasks))

