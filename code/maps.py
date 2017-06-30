import matplotlib.colors as colors
import matplotlib.cm as cmx
from mpl_toolkits.basemap import Basemap as Basemap
from matplotlib.colors import rgb2hex
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.patches import PathPatch
map_path = 'map_shapes/gz_2010_us_050_00_5m'
empire_file = 'map_shapes/empire_info.pkl'
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle

def load_maps():

    m = Basemap(llcrnrlon=-119,llcrnrlat=22,urcrnrlon=-64,urcrnrlat=49,
                projection='lcc',lat_1=33,lat_2=45,lon_0=-95,resolution='h')

    #m.drawcountries(linewidth=0.5, linestyle='solid', color='w',antialiased=1, ax=None, zorder=None)

    shp_info = m.readshapefile(map_path,'counties',drawbounds=False)
    emp = pickle.load(open(empire_file))
    return m, emp #
    

def mplmap(m, empire_list, value, symmetric=True,topq=100, state_level=False,rounddig=0,fa = (None,None), plot_empire=True,no_colorbar = False, crange=None ):
    (empire, rectE, rectW, rectN, rectS) = empire_list
    value = value[~pd.isnull(value)]

    cmap = plt.get_cmap('coolwarm')
    #v = min(np.abs(value.min()), value.max())
    v = min(np.abs(np.percentile(value, 100-topq)), np.percentile(value,topq))
    cNorm  = colors.Normalize(vmin=-1*v, vmax=v)
    crange = (-1*v, v)
    if not symmetric:
        cmap = plt.get_cmap('copper')
        cNorm  = colors.Normalize(vmin=value.min(), vmax=np.percentile(value,topq))
        crange = (value.min(), np.percentile(value,topq))
        
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
    return mplmap2(m, empire_list, value,scalarMap, state_level=state_level, crange=crange,
                       rounddig=rounddig,fa=fa, plot_empire=plot_empire,no_colorbar=no_colorbar)


def mplmap2(m, empire_list, value, scalarMap, state_level=False, crange=None,rounddig=0 ,
                fa = (None,None), plot_empire=True,no_colorbar = False):
    f, ax = fa
    if f is None:
        f, ax = plt.subplots(1,figsize=(12,6)) #gca() # get current axes instance
    ncou = 0
    patches = []
    (empire, rectE, rectW, rectN, rectS) = empire_list
    value = value[~pd.isnull(value)]
    for info, shape in zip(m.counties_info, m.counties):

        fip = info['STATE'] + info['COUNTY']

        sta = int(info['STATE'])
        cou = int(info['COUNTY'])
        cname = '%02d%03d' % (sta,cou)
        if fip in empire:
            if plot_empire == False:
                continue

            shape = empire[fip]
        #color = #rgb2hex(colors[cname])
        lookup_by = fip
        if state_level: lookup_by = info['STATE']        
        if not lookup_by in value.index: continue        

        poly = Polygon(shape,facecolor=scalarMap.to_rgba(value.loc[lookup_by]),edgecolor='w',lw=0.05)
            #patches.append(poly)
        ax.add_patch(poly)
        
    import matplotlib.patches as patches
    from matplotlib.path import Path
    '''    
    codes = [Path.MOVETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.CLOSEPOLY,
             ]
        
    path = Path([(rectE,rectS), (rectE,rectN), (rectW, rectN), (rectW, rectS),(0,0)], codes)
    '''
    if plot_empire:

        codes = [Path.MOVETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.CLOSEPOLY,
                 ]

        mid = np.mean([rectN, rectS])
        path = Path([(rectE,mid),(rectE,rectS), (rectW,rectS), (rectW, mid),
                    (rectE, mid),(rectE,rectN),(rectW,rectN),(rectW,mid),(0,0)], codes)

        ax.add_patch(
            patches.PathPatch(path,
                              edgecolor='black', facecolor='none')
        )        
    # cycle through state names, color each one.
    ax.relim()
    ax.set_axis_off()
    # update ax.viewLim using the new dataLim
    ax.autoscale_view()
    scalarMap._A = []
    if not crange:
        crange = [value.min(), value.max()]
    
    cbticks = [round(i,rounddig) for i in np.linspace(crange[0],crange[1],5)]
    bbox = ax.get_position()
    axcb, cbticks = [None, None]
    if not no_colorbar:
        axcb = f.add_axes([bbox.xmax,bbox.ymin + .25*(bbox.ymax-bbox.ymin),.01*bbox.width,.25*bbox.height])
    	cbticks = [round(i,rounddig) for i in np.linspace(crange[0],crange[1],5)]
    
        cb = f.colorbar(scalarMap,  #shrink=.25,
                    ticks=cbticks, cax= axcb)
    
    return f,ax, axcb, cbticks
