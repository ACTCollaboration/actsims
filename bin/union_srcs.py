from __future__ import print_function
from orphics import maps,io,cosmology,stats
from pixell import enmap,reproject,utils
import numpy as np
import os,sys
from soapack import interfaces as sints
from scipy import spatial
import matplotlib.pyplot as plt

#version = "06192019"
version = "20190819"

def merge_duplicates(ras,decs, rlim=1*utils.arcmin):
    """Modified from Sigurd's enlib.dory. Given a point source catalog 
    which might contain duplicates, detect these duplicates
    and merge them to produce a single catalog with no duplicates. Sources are considered
    duplicates if they are within rlim of each other. Merging uses averaging always. 
    rlim should be adjusted
    to fit the exerpiment beam. The default is appropriate for ACT."""
    # Normalize positions first. This could miss some mergers on the edge.
    ras = utils.rewind(ras, 0)
    pos    = np.array([ras*np.cos(decs),decs]).T
    tree   = spatial.cKDTree(pos)
    groups = tree.query_ball_tree(tree, rlim)
    done   = np.zeros(len(ras),bool)
    ocat   = []
    for gi, group in enumerate(groups):
        # Remove everything that's done
        group = np.array(group)
        group = group[~done[group]]
        if len(group) == 0: continue
        # Nothing to do for groups with only one member
        if len(group) == 1:
            done[group[0]] = True
            ocat.append((ras[group[0]],decs[group[0]]))
        else:
            gras  = ras[group]
            gdecs  = decs[group]
            era = gras.mean()
            edec = gdecs.mean()
            ocat.append((era,edec))
            done[group] = True
    return np.array(ocat)


eras,edecs,snrs,sizes = sints.get_act_mr3f_extended_sources(sn_cut=90.,size_cut=0.0)
print("Number of bright extended sources : ", len(eras))
bras,bdecs = sints.get_act_mr3f_cut_sources()
print("Number of bright cut sources : ", len(bras))

debug = False
rlim = 1.0
jras = np.append(eras , bras)
jdecs = np.append(edecs , bdecs)
ocat = merge_duplicates(jras*utils.degree,jdecs*utils.degree, rlim=rlim*utils.arcmin)  / utils.degree

extras = [(185.0433, 2.0733),
          (181.922, -1.108),
          (166.408, 2.043)]

print(ocat.shape)
for extra in extras:
    ocat = np.vstack((ocat,extra))

print(ocat.shape)

io.save_cols("union_catalog_%s.csv" % version,(ocat[:,0],ocat[:,1]),delimiter=',',header='ra(deg),dec(deg) | Made using actsims/bin/union_srcs.py.',fmt='%.5f')
#ras,decs = sints.get_act_mr3f_union_sources()
#assert np.all(np.isclose(ras,ocat[:,0]))
#assert np.all(np.isclose(decs,ocat[:,1]))


if debug:
    xmin = -190 ;  xmax = 190
    ymin = -70 ; ymax = 35
    s = 0.1
    dpi = 100

    plt.scatter(jras,jdecs,s=s,marker=',')
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.savefig("scatter.png",dpi=dpi)
    plt.clf()

    for i,rlim in enumerate([0.000001,0.001,0.1,1,2,4,6,10]):
        ocat = merge_duplicates(jras*utils.degree,jdecs*utils.degree, rlim=rlim*utils.arcmin)  / utils.degree
        print(rlim,ocat.shape)
        plt.scatter(ocat[:,0],ocat[:,1],s=s,marker=',')
        plt.xlim(xmin,xmax)
        plt.ylim(ymin,ymax)
        plt.savefig("scatter_%d.png" % i,dpi=dpi)
        plt.clf()

