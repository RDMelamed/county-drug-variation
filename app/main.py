import numpy as np
import pandas.io.sql as psql
import sqlite3 as sql
from bokeh.layouts import layout, widgetbox
from bokeh.layouts import row,widgetbox, column
from bokeh.plotting import figure, output_file, show #, cursession
from bokeh.models.widgets import Slider, Select, TextInput
from bokeh.io import curdoc
from bokeh.sampledata.movies_data import movie_path
from bokeh.models import Range1d
from bokeh.models import ColumnDataSource, LinearColorMapper
from bokeh.models import HoverTool
from bokeh.io import output_notebook
from bokeh.io import show
from bokeh.palettes import brewer
import sys
sys.path.append('../code')
import clean

from itertools import chain
import sys
#sys.path.append('bokeh/bokeh')
from bokeh.io import show
from bokeh.models import (
    ColumnDataSource,
    HoverTool,
    LinearColorMapper,
    Div
)
import numpy as np
import pandas as pd
import copy
from os.path import dirname, join
import pickle
import components
print 'step ONE'
dat = pd.read_pickle('data/pca_results_bokeh.pkl') #svd_dat.pkl') #cc100xscores.pkl') #county_cat.pkl')

outers = pickle.load(open('data/outer_counties.pkl'))
from bokeh.sampledata.us_counties import data as counties
natct =  pd.read_pickle('data/natct.pkl')
    #show(pusa)
    
#from bokeh.sampledata.unemployment import data as unemployment
#ptrend.name('trend')
#pusa, hover = create_figurecounties(.3)
#pusa.name('usamap')

dimlist = list(dat['xscores'].columns)
curdim = dimlist[0]
value = dat['xscores'][curdim]
VALLIST = list(value)
ix = [(int(i[:2]),int(i[2:])) for i in list(value.index)]
countysource = ColumnDataSource(data = {'x':[counties[k]['lons'] if not k[0] in [2,15] else outers[k]['lons'] for k in ix],
                    'y':[counties[k]['lats'] if not k[0] in [2,15] else outers[k]['lats'] for k in ix],
                    'value':VALLIST, #list(value),
                    'cty':list(value.index),
                   'name':list(natct.loc[value.index,'name'])})
pusa, colormapper = components.plot_usa_value(countysource)
pusa.title.text =value.name + ' Projection'

def update(addr, olddim, new):
    if olddim == new: return
    #dim = dimselect.value
    print 'changing to ' + new
    global curdim
    curdim = new
    global colormapper

    colormapper = components.get_cmap(dat['xscores'][curdim])
    #global pusa
    #global countysourceData
    #countysourceData = countysource.data.copy()
    #countysourceData['value'] = dat['xscores'][curdim]

    print 'yep... changed to ' + curdim
    print '--'.join(countysource.data.keys())
    global VALLIST
    VALLIST = list(dat['xscores'][new])
    countysource.data['value'] = VALLIST
    
    #countysource.data['value'] = list(dat['xscores'][curdim])
    #= {k:(v if not k=='value' else list(dat['xscores'][curdim]))
    #                            for k,v in countysource.data.items()} #countysourceData

    print 'for ' + new + ' SHOULD BE now:' , countysource.data['value'][:5]
    print countysource.data['value'][:5]
    pusa.title.text = curdim + ' projection'
    #return

    print 'changed title.. ' + new 
    #pdat, ldat = get_cor_datas(new, dat)
    #corpoints.data = pdat
    #corlines.data = ldat
    corfig.yaxis.axis_label = 'correlation with ' + new
    

    points, lines, demtexts, drugtexts, ymax, demopt, drugopt = components.get_cor_datas(new, dat)
    print 'got new cordata', (1.2*ymax[0], 1.2*ymax[1])
    #global yrange
    #yrange = Range1d(1.2*ymax[0], 1.2*ymax[1])
    #yrange = Range1d(-2,2)
    #cursession().store_objects(p) 
    corfig.y_range.start = 1.2*ymax[0]
    corfig.y_range.end = 1.2*ymax[1]     

    print 'assigning new cordata'
    cordata['points'].data = points
    cordata['lines'].data = lines
    
    cordata['demtext'].data = demtexts
    
    cordata['drugtext'].data = drugtexts
    print 'assigned texts , doing ymax', ymax[1]    
    #print 'ymax = ', ymax
    dcopy = cordata['xlabeldrug'].data.copy()
    dcopy['y'] = [ymax[1]]
    print dcopy
    cordata['xlabeldrug'].data = dcopy
    
    print 'assigned drug..', ymax[1]
    dcopy = cordata['xlabeldem'].data.copy()
    dcopy['y'] = [ymax[1]]
    cordata['xlabeldem'].data = dcopy
    #cordata['xlabeldrug'].data['y'] = ymax[1]
    print 'changing select...'
    drugselect.options = drugopt
    demselect.options = demopt
    global curdem
    curdem = demopt[0]
    curdrug = drugopt[0]
    drugselect.value = curdrug
    demselect.value = curdem    
    pairplots['dimdem'].xaxis.axis_label = new
    pairplots['dimdrug'].xaxis.axis_label = new    
    pairplots['dimdem'].yaxis.axis_label = clean.mystrip(curdem)
    pairplots['dimdrug'].yaxis.axis_label = clean.mystrip(curdrug)
    pairplots['demdrug'].yaxis.axis_label = clean.mystrip(curdrug)
    pairplots['demdrug'].xaxis.axis_label = clean.mystrip(curdem)
    logbutton.value='Linear scale'
    print 'done changing select...'
    #pair_update()    
    

def pair_update():
    r = components.pair_sources(dat['xscores'][curdim],demodat[curdem], drugdat[curdrug])
    for k, v in r.items():
        pairsources[k].data = v
    
        
        
def update_drug(addr, cur, new):
    if cur == new: return
    global curdrug
    curdrug = new
    pairplots['dimdrug'].yaxis.axis_label = clean.mystrip(new)
    pairplots['demdrug'].yaxis.axis_label = clean.mystrip(new)    
    pair_update()
    print 'here is a question:' +  curdim
    print countysource.data['value'][:5]

def update_dem(addr, cur, new):
    if cur == new: return
    global curdem
    curdem = new
    pairplots['dimdem'].yaxis.axis_label = clean.mystrip(new)    
    pairplots['demdrug'].xaxis.axis_label = clean.mystrip(new)    
    pair_update()    

from bokeh.models import BasicTicker, LogTicker
from bokeh.models import LinearAxis, LogAxis    
def update_scale(attr, old, new):
    ticker = LinearAxis()
    if new == "Log scale":
        ticker = LogAxis()
    pairplots['dimdem'].yaxis[0] = ticker
    #pairplots['demdrug'].xaxis[0].ticker = ticker    

desc = Div(text=open(join(dirname(__file__), "description.html")).read(), width=350)

dimselect = Select(title="Component", value=dimlist[0], options=dimlist)
dimselect.on_change('value', lambda attr, old, new: update(attr, old, new))
#countysourceData = countysource.data
corfig, cordata, dems, drugs = components.corrplot(dimlist[0], dat)
curdrug = drugs[0]
drugselect = Select(title="Compare drug", value=curdrug, options=drugs)
drugselect.on_change('value', lambda attr, old, new: update_drug(attr, old, new))
curdem = dems[0]
demselect = Select(title="Compare demographics", value=curdem, options=dems)
demselect.on_change('value', lambda attr, old, new: update_dem(attr, old, new))

curlogdem = False
logbutton = Select(title="", value="Linear scale", options=["Linear scale","Log scale"])
logbutton.on_change('value', lambda attr, old, new: update_scale(attr, old, new))

drugdat = pd.read_pickle('data/partition_county.pkl')
demodat = pd.read_pickle('data/demo_deredundanted.pkl')
pairplots, pairsources = components.create_pair_plots(dat['xscores'][curdim],
                                        demodat[curdem], drugdat[curdrug])
totlay = column(row(column(desc, dimselect), pusa),
                row(corfig, column(drugselect,
                                    demselect,
                                       #row(demselect, logbutton),
                                    pairplots['dimdem'], pairplots['dimdrug'],
                                    pairplots['demdrug'])))
#p, h = create_figure(.1)
curdoc().add_root(totlay)
curdoc().title = "CCA"
