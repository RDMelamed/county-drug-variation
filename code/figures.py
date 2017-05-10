import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from itertools import chain
import pdb
from scipy import stats

def bipartite_unlabeled(dim, cutdim, cutcor, cc,drugcolors,demcolors, c2d=None):
    f,ax = plt.subplots(1)
    f.set_size_inches(15,12)
    demcor = 'dem_v'
    topl = cc[demcor].loc[np.abs(cc[demcor][dim]) > cutdim, dim].sort_values()
    
    texts = []
    topdrchar = cc['drug_u'].loc[np.abs(cc['drug_u'][dim]) > cutdim, dim].sort_values()
    #topdrchar = topdrchar.drop('brand',axis=0)
    if not c2d is None:
        sel = np.tile(True, topdrchar.shape)
        for (i, c) in enumerate(topdrchar.index):
            if c in c2d:
                if len(c2d[c])==1: sel[i] = False
        print 'keeping', sel.sum()
        print 'rejecting', (~sel).sum()
        #print len(c2d['Bioflavonoids & Comb'])
        topdrchar = topdrchar[sel]
    x2 = 7
    ax.set_xlim([-.5,x2 + .5])
    #ax.set_ylim([-.85,.85])
    ax.set_xticks([1,x2])
    ax.set_xticklabels(['demographic\ncharacteristic','drug\ncharacteristic'])
    
    #ax.set_axis_off()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylabel('correlation with ' + dim)
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    x = cc['drug_dem'].loc[topdrchar.index, topl.index].stack()
    x = x[x > cutcor]

    for i in x.index:
        ax.plot([1.05,6.95],[topl.loc[i[1]],topdrchar.loc[i[0]]],color='c',
                linewidth=15*x.loc[i],alpha=.2)
    ax.plot(ax.get_xlim(),[0,0],'--',color=[.6,.6,.6])
    for i in topl.index:
        ax.plot(1, topl.loc[i],'.',markersize=15, color=demcolors[i],
                    markeredgewidth=.5, markeredgecolor=[.5,.5,.5])
    for i in topdrchar.index:
        ax.plot(x2, topdrchar.loc[i],'.',markersize=15,
                    color='k' if not i in drugcolors else drugcolors[i])
    return f,ax

def xy(x,y,ax,highlight=None):
    ax.plot(x['dat'],y['dat'],'.',color='gray',markersize=3,label='_nolegend_') #'Access to exercise opportunities Value'
    if highlight:
        for col, sel in highlight.items():
            ax.plot(x['dat'].loc[sel],y['dat'].loc[sel],'.',color=col,markersize=4)
    ax.set_xlabel(x['lab'])
    if 'log' in x:
        ax.set_xscale('log')
    else:
        ax.locator_params(axis='x',tight=True, nbins=5)
    if 'log' in y:
        ax.set_yscale('log')
    else:
        ax.locator_params(axis='y',tight=True, nbins=5)
    ax.set_ylabel(y['lab']) #'Repl Preps, Potassium Supp deviance')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    #ax.locator_params(nticks=4, axis='x')
    #ax.xaxis.major.locator.set_params(nticks=4) 
    #ax.locator_params(nticks=4) #,axis='x')
    ax.set_title('cor={:1.2f}'.format(stats.spearmanr(x['dat'],y['dat']).correlation))

def mystrip(i):
    return i.replace('Agents','').replace('Value','').replace('Percent of population','').replace('that is ','').replace('brandred','brand').replace('startdate_norm','release_date').replace('nadac','nadac ($)').strip().replace('otcred','OTC').replace('2010','').replace('2011','').replace('Female','F').replace('Male','M').replace('prevalence','').replace('(years)','').replace('(%)','').replace(', *','').replace('Native Hawaiian or Other Pacific Islander','Hawaiian/Pacific Islander').replace('Drugs','').replace('Drug','').replace('Misc','')

def circos_files(dim, cutdim, cutcor, cc,democlus, drugcolor,redname,bluename, nm=None, c2d=None):
    demcor = 'dem_v'
    topl = cc[demcor].loc[np.abs(cc[demcor][dim]) > cutdim, dim].sort_values()
    texts = []
    topdrchar = cc['drug_u'].loc[np.abs(cc['drug_u'][dim]) > cutdim, dim].sort_values()
    print 'dropping eyedrops'
    topdrchar = topdrchar.drop('Eyewash/Eyestrm/Lubr/Tear',axis=0)    
    highlights = open('circos_files/highlights.txt','w')
    for demval in topl.index:
        segment = 'demblue' if topl.loc[demval] < 0 else 'demred'
        val = np.abs(int(topl.loc[demval]*100))
        highlights.write('\t'.join([segment, str(val), str(val+1), mystrip(demval),
                                       'fill_color=' +democlus[demval] + ',' + 
                                      'color='+democlus[demval]]) + '\n')
    for drugval in topdrchar.index:
        segment = 'drugblue' if topdrchar.loc[drugval] < 0 else 'drugred'
        val = np.abs(int(topdrchar.loc[drugval]*100))
        colors = ''
        if drugval in drugcolor:
            colors = 'fill_color=' + drugcolor[drugval] + ',color=' + drugcolor[drugval]
        highlights.write('\t'.join([segment, str(val), str(val+1),
                                       mystrip(drugval.replace('Antitussives/Cold Comb','Allergy/Cold'))]) + '\t' + colors + '\n')
    x = cc['drug_dem'].loc[topdrchar.index, topl.index].stack()
    x = x[x > cutcor]
    links = open('circos_files/links.txt','w')
    for i in x.index:
        demval = int(topl.loc[i[1]]*100)
        drugval = int(topdrchar.loc[i[0]]*100)
        demsegment = 'demblue' if demval < 0 else 'demred'
        drugsegment = 'drugblue' if drugval < 0 else 'drugred'
        demval = np.abs(demval)
        drugval = np.abs(drugval)
        colors = ''
        if i[0] in drugcolor: colors = ',color=' + drugcolor[i[0]]
        links.write('\t'.join([demsegment, str(demval), str(demval + 1),
                            drugsegment, str(drugval), str(drugval + 1)])
                        + '\t' + 'thickness=' + str(int(20*x.loc[i])) + 'p' + colors + '\n')
    highlights.close()
    links.close()
    segments = open('circos_files/segments.txt','w')
    segments.write('\t'.join(['chr','-','drugred','drugs-' + redname,
                                  str(int(cutdim*100)-1), str(int(topdrchar.max()*100) + 1),'vlred'] )+ '\n')
    segments.write('\t'.join(['chr','-','drugblue','drugs-' + bluename,
                                  str(int(cutdim*100)-1), str(-1*int(topdrchar.min()*100) + 1),'vlblue']) + '\n')
    segments.write('\t'.join(['chr','-','demred','demographics-' + redname,
                                  str(int(cutdim*100)-1), str(int(topl.max()*100) + 1),'vlred']) + '\n')
    segments.write('\t'.join(['chr','-','demblue','demographics-' + bluename,
                                  str(int(cutdim*100)-1), str(-1*int(topl.min()*100) + 1) ,'vlblue']) + '\n')
    segments.close()
    
def img_only(axmatrix, fig, class_cor, colormap_string, cmin, cmax,
                xticks=True):
    cax = axmatrix.imshow(class_cor, #yidx
                          cmap=colormap_string,vmin=cmin, vmax = cmax, #-1*absmax, vmax=absmax,
                          interpolation='nearest',
                            aspect = 'auto') #,extent=[class_cor.shape[1]-1,0, class_cor.shape[1]-1,0]) #auto
    if xticks:
        axmatrix.set_xticks(range(class_cor.shape[1]))
        axmatrix.set_xticklabels(class_cor.columns,rotation=90)
    else:
        axmatrix.set_xticks([])
    axmatrix.set_yticks(np.array(range(class_cor.shape[0]))) #- .5)
    axmatrix.set_yticklabels(class_cor.index) #.iloc[yidx,:]
    

    axcbar = fig.add_axes([-.1, #0.95 if not no_row_dendrogram else .85 ,
                            0.1 if xticks else .1,0.08,0.3])
    axcbar.set_xticks([])
    axcbar.set_yticks([])
    #axcbar.set_axes([])
    axcbar.set_axis_off()
    cbar = fig.colorbar(cax, ax=axcbar) #rientation='horizontal',def

def mycont(x, y, scores, xryr, rl, bl):
    perc = np.percentile(scores,[10,90])
    sel = (scores < perc[0]) | (scores > perc[1])
    f, ax = plt.subplots(1)
    v = min(np.abs(np.percentile(scores, 0)), np.percentile(scores,100))
    crange = (-1*v, v)
    print crange
    sc = ax.scatter(x.loc[sel], y.loc[sel],s=4, cmap='bwr',vmin=crange[0],vmax=crange[1],
                c = scores.loc[sel],linewidth=.1)
    f.colorbar(sc,shrink=.4,ticks=[np.round(i,1) for i in np.linspace(scores.min(),scores.max(),5)])
    X, Y = np.mgrid[xryr[0]:xryr[1]:100j,xryr[2]:xryr[3]:100j,]
    positions = np.vstack([ X.ravel(), Y.ravel()])
    #return
    reds = scores.loc[scores > perc[1]].index
    kernel = stats.gaussian_kde(pd.concat((x.loc[reds], y.loc[reds]),axis=1).transpose())
    Z = np.reshape(kernel(positions).T, X.shape)
    ax.contourf(X, Y, Z, 6, cmap='Reds', levels=rl, alpha = .75) #,linewidths=3)

    blues = scores.loc[scores < perc[0]].index
    kernel = stats.gaussian_kde(pd.concat((x.loc[blues], y.loc[blues]),axis=1).transpose())
    Z = np.reshape(kernel(positions).T, X.shape)
    ax.contourf(X, Y, Z, cmap='Blues', levels=bl, alpha = .6) #extend='neither') #,linewidths=3)
    ax.set_xticks(np.linspace(-1,1,3))
    ax.set_yticks(np.linspace(-3,3,3))

    f.set_size_inches(3,2.5)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    return f, ax    
    
def drug_trend_images(modc, drug, todo):
    f = plt.figure(figsize=(8,3))
    axen = [(f.add_axes([0,.1,0.4,.9]), f.add_axes([.35,.2,0.1,.6])),
            (f.add_axes([.6,.1,0.4,.9]),f.add_axes([.95,.2,0.1,.6]))] 
    #axen = plt.subplots(1,len(todo))
    for axpair, yv in zip(axen, todo):
        #
        ax, axcbar = axpair
        regroup = todo[yv]
        drugdat = modc.transpose().stack('type').transpose()[drug].reset_index()
        #sel = np.ones(drugdat.shape[0], dtype=bool)
        sel = (drugdat['age'] < 70) & (drugdat['year'] < 11) & (drugdat['visits'] < 15)
        for dim in regroup:
            sel = sel & (drugdat[dim] == regroup[dim])
        #print sum(sel)
        drugdat = drugdat.loc[sel,:]    
        #set(drugdat.columns)-set(['count','denom']) - set( + todo.keys() )
        drugdat = drugdat.groupby(['year',yv]).sum().reset_index()
        drugdat['rate'] = drugdat['count']/drugdat['denom'].map(float)        
        #drugdat = drugdat.loc[~np.isnan(drugdat['rate']),:]
        drugdat.loc[np.isnan(drugdat['rate']),'rate'] = 0
        ratmat = drugdat.loc[:,['year',yv,'rate']].set_index(['year',yv]).transpose()
        #.stack('age')
        #print 'here', ratmat.shape    
        #print ratmat.iloc[:5,:3]

        #print ratmat.iloc[:5,:].index

        ratmat = ratmat.stack(yv)
        forimg = ratmat
        vmax = np.percentile(forimg,99)
        cax = ax.imshow(forimg,cmap='hot',interpolation='nearest', aspect = 'auto',vmax=vmax)
        ax.set_xlabel('year',fontsize=14)
        ax.set_ylabel(yv.replace('visits','number Rx'),fontsize=14)
        ages = np.array(forimg.index.levels[1])
        atick = range(0,len(ages),2)
        ax.set_yticks(atick)
        tit = ';'.join([dim.replace('visits','number Rx') + '=' + str(val) for dim, val in regroup.items()
                                         if not dim=='nbin'])
        if len(tit) > 0:
            ax.set_title('among: ' + tit)
        ax.set_yticklabels(ages[atick])
        yrs = np.array(forimg.columns)
        yrtick = range(0,len(yrs),2)
        ax.set_xticks(yrtick)
        ax.set_xticklabels(yrs[yrtick] + 2003, rotation = 90)
        v90 = np.percentile(forimg,90)
        v10 = np.percentile(forimg,10)
        ticky=[np.round(v10,3),np.round(np.mean(forimg.mean()),3),np.round(v90,3)]
        #print ticky
        axcbar.set_xticks([])
        axcbar.set_yticks([])
        #axcbar.set_axes([])
        axcbar.set_axis_off()
        cbar = f.colorbar(cax, ax=axcbar, ticks=ticky) #rientation='horizontal',
    f.text(.5,1.1,drug + ' / person-year',ha='center',va='top',fontsize=14)
    return f

def make_compress(unstacked, cuts, dimname):
    groupers = list(set(['year','age','nbin','visits','type']) - set([dimname]))
    tocat = []
    unstacked = unstacked.reset_index()                 
    for ix in range(len(cuts)):
        sel = unstacked[dimname] >= cuts[ix]
        if ix < len(cuts)-1: sel = sel & (unstacked[dimname] < cuts[ix + 1])
        sub = unstacked.loc[sel,:].groupby(groupers).sum()
        sub[dimname] = cuts[ix]
        sub.set_index(dimname,append=True,inplace=True)
        tocat.append(sub)
    unstacked = pd.concat(tocat, axis=0)
    unstacked = unstacked.reorder_levels(['year','age','visits','nbin','type'])
    unstacked = unstacked.sort_index() 
    return unstacked
