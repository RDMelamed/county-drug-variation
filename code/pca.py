import pandas as pd
import numpy as np
from itertools import chain
import pdb
from scipy import stats
from sklearn.utils.extmath import randomized_svd
def do_svd(devin, demo0, fda_nadac, c2d):
    ### Run PCA
    dev = (devin - devin.mean(axis=0))/devin.std(axis=0)
    U, Sigma, VT = randomized_svd(dev.transpose(),
                                      n_components=8,
                                            n_iter=7,
                                            random_state=42)
    # get projection of counties
    proj = pd.DataFrame((VT.transpose()*Sigma),index=dev.index)
    proj = proj.rename(columns={i:'V' + str(i+1) for i in proj.columns})
    ret = {'xscores':proj}
    ### Compile correlations

    # correlation of drug deviance with projection
    totcol = list(proj.columns) + list(dev.columns)
    ret['drugAll_u'] = pd.DataFrame(stats.spearmanr(proj,dev).correlation,
                      index = totcol, columns = totcol).loc[dev.columns, proj.columns]
    
    # correlation of demographics with drug classes / characteristics
    charsel = ['acute','chronic','brandred','startdate_norm','nadac']    
    demParti = pd.DataFrame(0,index=c2d.keys() + charsel,columns=demo0.columns)
    drugu = pd.DataFrame(0, index =c2d.keys() + charsel, columns = proj.columns)
    for c in c2d.keys():
        partitionAcrossCounties = dev.loc[:,dev.columns.isin(c2d[c])].median(axis=1)
        demParti.loc[c,:] = stats.spearmanr(partitionAcrossCounties, demo0).correlation[0,1:]
        drugu.loc[c,:] = stats.spearmanr(partitionAcrossCounties, proj).correlation[0,1:]

    # silly correlation of correlations... probably useless
    for c in charsel:
        sel = ~pd.isnull(fda_nadac[c]) & (fda_nadac.index.isin(dev.columns))
        v = fda_nadac.loc[sel,c]
        county_char = stats.spearmanr(v, dev.loc[:, v.index].transpose()).correlation[0,1:]
        demParti.loc[c,:] = stats.spearmanr(county_char, demo0).correlation[0,1:]
        drugu.loc[c,:] = stats.spearmanr(v,ret['drugAll_u'].loc[v.index, :]).correlation[0,1:]

    ret['drug_dem'] = demParti        
    ret['drug_u'] = drugu
    
    # correlation of demographics with projections
    totcol = list(proj.columns) + list(demo0.columns)
    ret['dem_v'] = pd.DataFrame(stats.spearmanr(proj,demo0).correlation,
                      index = totcol, columns = totcol).loc[demo0.columns, proj.columns]
    return ret
'''
    svdhypbeta = cca_viz.ccaize4(dev50.loc[demo0.index,drugset], 
                             demo0, svd50,fda_nadac, charDem.loc[:,charsel],
                                [c for c,v in c2d.items() if len(set(v) & set(drugset))>0], c2d)
    
    return proj, U, Sigma, VT
'''
