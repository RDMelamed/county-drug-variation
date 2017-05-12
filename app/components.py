from bokeh.models.glyphs import MultiLine, Circle, Text
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, LinearColorMapper
from bokeh.models import HoverTool
from bokeh.palettes import brewer
from bokeh.models import FixedTicker
from bokeh.core.enums import TextAlign

import numpy as np
import pandas as pd
import sys
sys.path.append('../code')
import clean
def get_cmap(value):
    #v = max(np.abs(value.min()), np.abs(value.max()))
    topq = 99
    v = min(np.abs(np.percentile(value, 100-topq)), np.percentile(value,topq))
    return LinearColorMapper(palette=brewer['RdYlGn'][11],
                                                        low = -1*v,
                                                        high= v)
    
def plot_usa_value(source):
    '''
    ix = [(int(i[:2]),int(i[2:])) for i in list(value.index)]
    countysource = {'x':[counties[k]['lons'] if not k[0] in [2,15] else outers[k]['lons'] for k in ix],
                    'y':[counties[k]['lats'] if not k[0] in [2,15] else outers[k]['lats'] for k in ix],
                    'value':list(value),
                    'cty':list(value.index),
                   'name':list(natct.loc[value.index,'name'])}

    source = ColumnDataSource(data=countysource)
    '''
    pusa = figure(tools="hover,box_zoom,save,reset", x_axis_location=None, y_axis_location=None)

    colormapper = get_cmap(source.data['value'])
    #return pusa, source, colormapper    
    pusa.patches('x', 'y', source=source,
              fill_color={'field':'value',
                          'transform':colormapper},
              line_color="white", line_width=0.05) #,alpha={'field':'missing'}

    pusa.grid.grid_line_color = None
    hover = pusa.select_one(HoverTool)
    hover.point_policy = "follow_mouse"
    hover.tooltips = [
        ("value", "@value"),
        ("cty", "@cty"),
        ("name","@name")
    ]
    #pusa.title.text =value.name + ' Projection' #'{:s}\tyear={:d}\tage={:d}\tyear_Rx={:d}\tyears_obs={:d}\n'.format(drug,year,age,vis,nyear)    
    return pusa, colormapper

def get_texts(topl, spac, xval):
    curset = []
    prvy = -2
    alt = True
    ys = []
    texts = []
    for i in topl.index:
        #texts.append(ax.text(1, topl[i], i))
        
        if topl[i] < prvy + spac and len(' | '.join(curset)) <  60:
            curset.append(i)
        else:
            if len(curset)> 0:
                #if ix < len(curset): tojoin.append(curset[ix]) # = curset
                #lab = clean.mystrip('\n'.join(tojoin))
                ys.append(topl.loc[curset].mean())
                texts.append(' | '.join([clean.mystrip(v) for v in curset]))
            curset = [i] 
        prvy = topl[i]
    texts.append(' | '.join([clean.mystrip(i) for i in curset]))
    ys.append(topl.loc[curset].mean())
    #return ys, texts

    #def text_dict(xval, topl):
    #ys, texts = get_texts(topl, .03)
    aligns = [TextAlign.left,TextAlign.right]*(len(texts)/2)
    aligns = aligns[:len(ys)]
    xs = [xval]*len(ys)
    return {'y':ys, 'x':xs,'text':texts}
x1 = 1 #'demographic'
x2 = 2 #'drug'
from bokeh.models.formatters import TickFormatter

                
def get_cor_datas(dim, cc, cutdim=.18, cutcor=.2):
    topl = cc['dem_v'][dim]
    topl = topl.loc[np.abs(topl) > cutdim].sort_values()
    topdrchar = cc['drug_u'][dim]
    topdrchar = topdrchar.loc[np.abs(topdrchar) > cutdim].sort_values()


    points = {'xs':[1]*len(topl) + [x2]*len(topdrchar),
                    'ys':list(topl) + list(topdrchar),
                    'name':list(topl.index) + list(topdrchar.index),
                    'correlation':list(topl) + list(topdrchar)}
    x = cc['drug_dem'].loc[topdrchar.index, topl.index].stack()
    x = x[x > cutcor]
    
    lines = {'xs':[[x1,x2]]*len(x),
            'ys':[[topl.loc[i[1]],topdrchar.loc[i[0]]] for i in x.index],
            #'name':[i[1] + ' - ' + i[0] for i in x.index],
            #'correlation':list(x),
                 'ws':list(9*x)}
    spac = .03
    demdict = get_texts(topl, spac, x1)
    drugdict = get_texts(topdrchar, spac, x2)
    #for k, v in drugdict.items():
    #    demdict[k].append(v)
        
    return (points, lines, demdict, drugdict, [min(topl.min(), topdrchar.min()),
                                                   max(topl.max(), topdrchar.max())],
                list(topl.index)[::-1], list(topdrchar.index)[::-1])


def corrplot(dim, cc, cutdim=.18, cutcor=.2):
    points, lines, demtexts, drugtexts, ymax, dems, drugs = get_cor_datas(dim, cc, cutdim, cutcor)
    linesource = ColumnDataSource(data=lines)
    pointsource = ColumnDataSource(data=points)    
    demtext = ColumnDataSource(data=demtexts)
    drugtext = ColumnDataSource(data=drugtexts)        
    f = figure(width=500,plot_height=620,
                   tools='hover,box_zoom,save,reset,pan',x_axis_location=None,
                   x_range=(x1-5, x2+5), y_range=(1.2*ymax[0], 1.2*ymax[1]))
    f.yaxis.axis_label = 'correlation with ' + dim
    f.grid.grid_line_color = None
    f.toolbar.logo = None
    gl = MultiLine(xs="xs",ys="ys",line_color='#00ffff',line_width='ws') #, source=linesource)
    f.add_glyph(linesource, gl)
    circles = f.circle(x="xs",y="ys",color='#003300',size=5, source=pointsource)
    #f.add_glyph(pointsource, circles)

    #print 'yrange max', f.y_range[1]
    dmgl = Text(x="x",y="y", text="text", text_align='right',x_offset=-2,text_font_size='8pt')
    drgl = Text(x="x",y="y", text="text", text_align='left',x_offset=2,text_font_size='8pt')
    '''
    for i, v in demtexts.items():
        print i, len(v)
    for i, v in drugtexts.items():
        print i, len(v)
    '''
    f.add_glyph(demtext, dmgl)
    f.add_glyph(drugtext, drgl)        
    print 'ymax',ymax
    demlabelsource = ColumnDataSource(data = {'x':[x1], 'y':[ymax[1]],
                                               'text':['Demographics']})
    druglabelsource = ColumnDataSource(data = {'x':[x2], 'y':[ymax[1]],
                                               'text':['Drugs']})
    demlabel = Text(x="x",y="y", text="text",
                      text_align='right',y_offset=-15,
                        text_font_size='10pt',text_font_style='bold')
    f.add_glyph(demlabelsource, demlabel)
    druglabel = Text(x="x",y="y", text="text",
                      text_align='left',y_offset=-15,text_font_size='10pt',
                         text_font_style='bold')
    f.add_glyph(druglabelsource, druglabel)                
    hover = f.select_one(HoverTool)
    hover.point_policy = "snap_to_data"
    hover.tooltips = [("name", "@name"),
                          ("correlation","@correlation")]
    hover.renderers = [circles]

    f.yaxis.minor_tick_line_color = None
    #f.xaxis.ticker = FixedTicker(ticks=[, x2])
    #f.xaxis.formatter = CategoricalTickFormatter(format="0.0%")
    return f, {'points':pointsource, 'lines':linesource,
                   'demtext':demtext, 'drugtext':drugtext, 'xlabeldem':demlabelsource,
                   'xlabeldrug':druglabelsource}, dems, drugs

def create_pairfig(sourcedat, xname, yname):
    f = figure(width=150,plot_height=150,
                   tools='')
    f.yaxis.minor_tick_line_color = None
    f.xaxis.minor_tick_line_color = None        
    f.yaxis.axis_label = clean.mystrip(yname)
    f.xaxis.axis_label = clean.mystrip(xname)
    f.xaxis.axis_label_text_font_size = "8pt"
    f.yaxis.axis_label_text_font_size = "8pt"    
    f.grid.grid_line_color = None
    f.toolbar.logo = None
    f.toolbar_location = None
    f.grid[0].ticker.desired_num_ticks = 3
    f.grid[1].ticker.desired_num_ticks = 5
    #f.yaxis.desired_num_ticks = 5   
    circles = f.circle(x="x",y="y",color='#003300',size=1, source=sourcedat)
    return f

def pair_sources(dimdat, demdat, drugdat):
    return {'dimdem':{'x':list(dimdat), 'y':list(demdat)},
            'dimdrug':{'x':list(dimdat), 'y':list(drugdat)},
            'demdrug':{'x':list(demdat), 'y':list(drugdat)}}
            
    
def create_pair_plots(dimdat, demdat, drugdat):


    dicts = pair_sources(dimdat, demdat, drugdat)
    sources = {k:ColumnDataSource(data=v) for k, v in dicts.items()}
    figures = {'dimdem':create_pairfig(sources['dimdem'],
                                           dimdat.name, demdat.name),
                'dimdrug':create_pairfig(sources['dimdrug'],
                                           dimdat.name, drugdat.name),
                'demdrug':create_pairfig(sources['demdrug'],
                                           demdat.name, drugdat.name)}
    return figures, sources
