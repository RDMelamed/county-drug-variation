import numpy as np
import pandas as pd
from itertools import chain
import pdb
from scipy import stats
    
def correlations(c2d, svddfilt, fda_nadac, cutndrug=8):
    pimp = {}
    cdo = c2d.keys()
    dim2 = 'brandred'
    for dim in svddfilt['xscores'].columns[:4]:
        
        xval = svddfilt['drugAll_u'][dim]
        y = fda_nadac.loc[xval.index,dim2] #noop['drugAll_u']        
        forhists = pd.DataFrame(float('nan'),index=cdo,
                                columns = ['corr' + dim2, 'med' + dim, 
                                            'med' + dim2,'std' + dim2,'std' + dim,
                                               'n','spearman' + dim2,'abs' + dim])
        
        for c in cdo: 
            sel = c2d[c] &set(xval.index)
            if len(sel) >= cutndrug:
                #cor = stats.spearmanr(x.loc[sel], y.loc[sel]).correlation
                cor = stats.pearsonr(xval.loc[sel], y.loc[sel])[0]
                cor2 = stats.spearmanr(xval.loc[sel], y.loc[sel]).correlation
                forhists.loc[c,:] = [cor, xval.loc[sel].median(), y.loc[sel].median(), y.loc[sel].std(),xval.loc[sel].std(),len(sel),cor2, xval.loc[sel].abs().max()]
        pimp[dim] = forhists.loc[~pd.isnull(forhists['corr' + dim2]),:]
    corbs = pd.concat({d:i['corrbrandred'] for d, i in pimp.items()},axis=1)
    sporbs = pd.concat({d:i['spearmanbrandred'] for d, i in pimp.items()},axis=1)
    spabs = sporbs.copy()
    for i in spabs.index:
        for j in spabs.columns:
             if np.abs(sporbs.loc[i,j]) > np.abs(corbs.loc[i,j]):
                spabs.loc[i,j] = corbs.loc[i,j]
    return spabs
