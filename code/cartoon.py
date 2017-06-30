from palettable.colorbrewer.sequential import Greens_9,Reds_9,YlGn_9,YlOrRd_9,RdPu_9,PuBuGn_9,Oranges_9
from palettable.colorbrewer.sequential import Purples_9,OrRd_9,Greys_9,Blues_9,BuPu_9,PuRd_9,PuBu_9,GnBu_9
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from itertools import chain
from scipy import stats
def comp_assign(scores, scperc):
    ctyass = pd.DataFrame(0, index=scores.index, columns=scores.columns[:4])
    for d in ctyass.columns:
        ctyass.loc[scores[d] < scperc[d]['perc'][0], d] = -1
        ctyass.loc[scores[d] > scperc[d]['perc'][1], d] = 1
    a = pd.DataFrame(stats.zscore(scores.iloc[:,:4],axis=0), index = scores.index) #.
    for c in ctyass.loc[ctyass.abs().sum(axis=1) > 0,:].index:
        z = np.array(a.loc[c,:].abs())
        ctyass.loc[c,z != z.max()] = 0
    return ctyass, a

def assign2color(ctyass, scores,scperc,NC):
    toplot = pd.DataFrame(float('nan'),index=ctyass.index,columns = ['cid'])
    color_i = 0

    for d in scperc.keys():
        incomp = ctyass[d]!=0
        dimval = scores[d]
        bins = np.percentile(dimval.loc[incomp],np.linspace(0,100,NC+1))[:-1]
        for b in range(len(bins)):
            bsel = incomp & (dimval >= bins[b])
            if b < len(bins) - 1:
                bsel = bsel & (dimval < bins[b+1])
            toplot.loc[bsel] = color_i
            print '{:s}, col={:d}, n={:d}'.format(d, color_i,bsel.sum())
            color_i += 1  
    return toplot

def get_color_map():
    s = 3
    e = 7
    cset = list(chain.from_iterable([[col.colors[s], col.colors[e]] 
            for col in [Greens_9, Purples_9, Blues_9, Reds_9]]))
    cset += [[256]*3]

    cset = np.array(cset)/256.0
    from matplotlib import colors as mcolors
    cmap = mcolors.LinearSegmentedColormap.from_list(name='clus', 
                                          colors=cset,N=len(cset))
    import matplotlib.cm as cmx
    cNorm = mcolors.Normalize(vmin=0, vmax=len(cset))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
    return cset, scalarMap

import matplotlib.gridspec as gridspec
import maps

def gridit_hviolins(scperc, scores, ctyass, ctcolor, scalarMap,cset, dsn2, m, emp,pset, partition2drug, state_level=False,percolor=4):
    ctcolor[pd.isnull(ctcolor)] = len(cset) - 1 #white
    allcty = [c['STATE'] +  c['COUNTY'] for c in m.counties_info]
    mi = set(allcty) - set(ctcolor.index)
    tocol = pd.concat((ctcolor, pd.DataFrame(len(cset) - 1, index=mi,columns=['cid'])),axis=0)

    
    gs = gridspec.GridSpec(5, len(scperc))
    f = plt.figure(figsize=(10,10))
    xl = [float('inf'),-1*float('inf')]
    yl = [float('inf'),-1*float('inf')]
    axen = []
    for i, d in  enumerate(scperc.keys()):
        axi = plt.subplot(gs[3,i])
        tomap = scores[d].copy() # pd.Series(0,index=.index)
        tomap.loc[ctyass[d]==0] = 0 # make it white float('nan')
        #print len(mi)
        tomap = pd.concat((tomap, pd.Series(0, index=mi)),axis=0)
        #print 'tomap:', tomap.shape
        #print sum(pd.isnull(tomap))
        maps.mplmap(m,emp, tomap, fa=(f,axi),plot_empire=False,no_colorbar=True,state_level=state_level)
        axi.set_title(scperc[d]['name'],fontsize=10)
        axen.append(axi)
        xl = [min(xl[0], axi.get_xlim()[0]), max(xl[1], axi.get_xlim()[1])]
        yl = [min(yl[0], axi.get_ylim()[0]), max(yl[1], axi.get_ylim()[1])]
    for ax in axen: 
        ax.set_xlim(xl)
        ax.set_ylim(yl)
    axi = plt.subplot(gs[:3,:])
    fo = maps.mplmap2(m, emp, tocol['cid'],scalarMap,fa=(f,axi), plot_empire=False,
                      no_colorbar=True,state_level=state_level)
    bbox = axi.get_position()
    axcb = f.add_axes([bbox.xmax-.15*(bbox.xmax-bbox.xmin),bbox.ymin + .3*(bbox.ymax-bbox.ymin),.01*bbox.width,.25*bbox.height])
    spa = percolor/2.0
    cb = f.colorbar(scalarMap,  #shrink=.25,
                    ticks=np.linspace(spa,len(cset)-1-spa,4), cax= axcb)
    axcb.set_yticklabels([scperc[d]['name'] for d in scperc])
    
    for i,d in enumerate(scperc.keys()):
        print i
        axi = plt.subplot(gs[4,i])
        if i == 0:
            axi.set_yticklabels(pset) #,rotation=45,ha='left')
        xval = range(len(pset))

        bs = axi.violinplot([list(dsn2.loc[scores[d] > scperc[d]['perc'][1],
                                     partition2drug[p]].median(axis=1)) for p in pset], 
                       positions = xval,showmeans=False,showextrema=False)
        for b in bs['bodies']: 
            b.set_color('r')
        bs = axi.violinplot([list(dsn2.loc[scores[d] < scperc[d]['perc'][0],
                                     partition2drug[p]].median(axis=1)) for p in pset], 
                       positions = xval,showmeans=False,showextrema=False)
        for b in bs['bodies']: 
            b.set_color('b')
        axi.set_xticks(xval)
        #if i > 1:
        axi.set_xticklabels(pset,rotation=90) #,rotation=45,ha='left')
        #else:
        #axi.set_xticklabels([]) #,rotation=45,ha='left')
        axi.spines['right'].set_visible(False)
        axi.spines['top'].set_visible(False)
        axi.set_yticks([-1,0,1])    
        axi.set_ylim(-2,2)          
    
    return f
