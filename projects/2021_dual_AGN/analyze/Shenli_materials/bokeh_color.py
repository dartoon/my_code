#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 15:59:38 2020

@author: tang
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import operator
#from sklearn.cluster import KMeans
from scipy.interpolate import interp1d

from bokeh.models import ColumnDataSource, Select, RangeSlider, Slider
from bokeh.layouts import row, column
from bokeh.plotting import figure
from bokeh.models import HoverTool,Range1d
from bokeh.io import output_file, show, curdoc

plt.rcParams['savefig.dpi'] = 200 
plt.rcParams['figure.dpi'] = 100 

data_path = '/Users/Dartoon/Downloads/shenli_data'
os.chdir(data_path)
sample = pd.read_csv('five_band_color.csv', index_col = 0)
sample['status'] = 'wait' # actually some has been observe, but not updated for the latest list yet
sample['telescope'] = 'wait' 
stellar = pd.read_csv('stellar.csv')
print(stellar)
#%%
stellar_color = np.array(stellar[['g-r','r-i']])
def distance(x,y):
    values = np.sqrt((x-stellar_color[:,0])**2+(y-stellar_color[:,1])**2)
    return min(enumerate(values), key=operator.itemgetter(1))
sample['qso_dis'] = sample.apply(lambda x: distance(x['qso_g-r'], x['qso_r-i'])[1], axis=1)
sample['com_dis'] = sample.apply(lambda x: distance(x['com_g-r'], x['com_r-i'])[1], axis=1)
sample['pos_qso_dis'] = sample.apply(lambda x: distance(x['qso_g-r'], x['qso_r-i'])[0], axis=1)
sample['pos_com_dis'] = sample.apply(lambda x: distance(x['com_g-r'], x['com_r-i'])[0], axis=1)

#%%
# Convert dataframe to column data source
qso_com = ColumnDataSource(data={
    'No' : sample.index,
    'z' : sample['Redshift'],
    'qso_y' : sample['qso_r-i'],
    'qso_x' : sample['qso_g-r'],
    'com_y' : sample['com_r-i'],
    'com_x' : sample['com_g-r'],
    'qso_dis' : sample['qso_dis'],
    'com_dis' : sample['com_dis'],
    'Sep' : sample['Sep(")'],
    'ID' : sample['ID'],
    'imgs' : sample['img'],
    'status' : sample['status']
})
#QSO = []
#done = []
#wait = []
#for index,item in enumerate(qso_com.data['status']):
#    if item == 'QSO':
#        QSO.append(index)
#    elif item == 'done':
#        done.append(index)
#    else:
#        wait.append(index)
#print(qso_com.data.index)
#print(qso_com.data[][QSO])
#%%
hover = HoverTool(
        tooltips="""
        <div>
            <div>
                <img
                    src="@imgs" height="150" width="150"
                    style="float: left; margin: 0px 15px 15px 0px;"
                    border="2"
                ></img>
            </div>
            <div>
                <span style="font-size: 17px; font-weight: bold;">@ID</span>
            </div>
            <div>
                <span style="font-size: 17px;">[@No]</span>
            </div>
        </div>
        """
    )

TOOLS = 'pan,box_zoom,box_select,lasso_select,reset'
    
# Create the figure: p
p1 = figure(y_axis_label='qso_r-i', x_axis_label='qso_g-r',tools=TOOLS)
#qso vs stellar locus
p1.circle('qso_x', 'qso_y', source=qso_com, color='grey', size=6, alpha = 0.3, hover_fill_alpha = 1.0, hover_fill_color = 'red')
p1.line(stellar['g-r'], stellar['r-i'], color = 'black', line_width = 1)
#p1.line((stellar['g-r'].iloc[sample['pos_qso_dis']], stellar['r-i'].iloc[sample['pos_qso_dis']]),('qso_x', 'qso_y'),source=qso_com,color = 'red',line_width = 1, alpha = 0, hover_fill_alpha = 1)

p2 = figure(y_axis_label='com_r-i', x_axis_label='com_g-r',tools=TOOLS)
p2.circle('com_x', 'com_y', source=qso_com, color='grey', size=6, alpha = 0.3, hover_fill_alpha = 1.0, hover_fill_color = 'blue')
p2.line(stellar['g-r'], stellar['r-i'], color = 'black', line_width = 1)

p3 = figure(y_axis_label='qso_dis', x_axis_label='qso_g-r',tools=TOOLS)
p3.circle('qso_x', 'qso_dis', source=qso_com, color='grey', size=6, alpha = 0.3, hover_fill_alpha = 1.0, hover_fill_color = 'red')

p4 = figure(y_axis_label='com_dis', x_axis_label='com_g-r',tools=TOOLS)
p4.circle('com_x', 'com_dis', source=qso_com, color='grey', size=6, alpha = 0.3, hover_fill_alpha = 1.0, hover_fill_color = 'blue')

p5 = figure(y_axis_label='qso_dis', x_axis_label='qso_r-i',tools=TOOLS)
p5.circle('qso_y', 'qso_dis', source=qso_com, color='grey', size=6, alpha = 0.3, hover_fill_alpha = 1.0, hover_fill_color = 'red')

p6 = figure(y_axis_label='com_dis', x_axis_label='com_r-i',tools=TOOLS)
p6.circle('com_y', 'com_dis', source=qso_com, color='grey', size=6, alpha = 0.3, hover_fill_alpha = 1.0, hover_fill_color = 'blue')

p2.x_range = p1.x_range = p4.x_range = p3.x_range = Range1d(-1,3)
p6.x_range = p5.x_range = Range1d(-1,3)
p2.y_range = p1.y_range = Range1d(-0.6,3)
p4.y_range = p3.y_range = Range1d(-0.2,2)
p6.y_range = p5.y_range = Range1d(-0.2,2)

plots = column(row(p1, p2),row(p3, p4),row(p5, p6))

p1.add_tools(hover)
p2.add_tools(hover)
p3.add_tools(hover)
p4.add_tools(hover)
p5.add_tools(hover)
p6.add_tools(hover)

range_slider = RangeSlider(start=0, end=24, value=(0,24), step=1, title='RA')
range_slider2 = RangeSlider(start=0, end=5, value=(0,5), step=0.1, title='Redshift')
range_slider3 = RangeSlider(start=0, end=4.2, value=(0,4.2), step=0.1, title='Sep(")')

# Create a dropdown Select widget: select
#options1 = ['total','close','medium','medium2','far']
#select1 = Select(title='Sep', options=options1, value='total')

options2 = ['all','wait','QSO','done','unclassified']
select2 = Select(title='status', options=options2, value='all')

# Define a callback function: update_plot
def update_plot(attr, old, new):
    # Read the current value of the slider: scale
    time_zone = range_slider.value
    redshift = range_slider2.value
    Sep = range_slider3.value
    status = select2.value
    # Update source with the new data values
    RA_con = (sample['RA'] >= time_zone[0]*15) & (sample['RA'] <= time_zone[1]*15)
    z_con = (sample['Redshift'] >= redshift[0]) & (sample['Redshift'] <= redshift[1])
    Sep_con = (sample['Sep(")'] >= Sep[0]) & (sample['Sep(")'] <= Sep[1])
    if status == 'all':
        sta_con = 1
    else:
        sta_con = (sample['status'] == status)
    com_con = RA_con & z_con & Sep_con & sta_con
    qso_com.data = {
                'No' : sample.loc[com_con].index,
                'z' : sample['Redshift'],
                'qso_y' : sample.loc[com_con, 'qso_r-i'],
                'qso_x' : sample.loc[com_con, 'qso_g-r'],
                'com_y' : sample.loc[com_con, 'com_r-i'],
                'com_x' : sample.loc[com_con, 'com_g-r'],
                'qso_dis' : sample.loc[com_con, 'qso_dis'],
                'com_dis' : sample.loc[com_con, 'com_dis'],
                'Sep' : sample.loc[com_con, 'Sep(")'],
                'ID' : sample.loc[com_con, 'ID'],
                'imgs' : sample.loc[com_con, 'img'],
                'status' : sample.loc[com_con, 'status']
                }

slider = Slider(start=10, end=50, value=10, step=1, title="bins")
qso_arr_hist, qso_edges = np.histogram(sample['qso_dis'], 
                           bins = 10, 
                           range = [0, 1.6])
com_arr_hist, com_edges = np.histogram(sample['com_dis'], 
                           bins = 10, 
                           range = [0, 1.6])
dis_hist_source = ColumnDataSource(data={
    'qso_num' : qso_arr_hist,
    'qso_left': qso_edges[:-1], 
    'qso_right': qso_edges[1:],
    'com_num' : com_arr_hist,
    'com_left': com_edges[:-1], 
    'com_right': com_edges[1:]})
def update_hist(attr, old, new):
    bins = slider.value
    qso_arr_hist, qso_edges = np.histogram(sample['qso_dis'], 
                               bins = bins, 
                               range = [0, 1.6])
    com_arr_hist, com_edges = np.histogram(sample['com_dis'], 
                               bins = bins, 
                               range = [0, 1.6])
    dis_hist_source.data = {
        'qso_num' : qso_arr_hist,
        'qso_left': qso_edges[:-1], 
        'qso_right': qso_edges[1:],
        'com_num' : com_arr_hist,
        'com_left': com_edges[:-1], 
        'com_right': com_edges[1:]
        }

p7 = figure(x_axis_label = 'qso_distance', y_axis_label = 'Numbers')
p7.quad(bottom=0, top='qso_num',left='qso_left', right='qso_right', source=dis_hist_source, fill_color='red', line_color='black')    
p8 = figure(x_axis_label = 'com_distance', y_axis_label = 'Numbers')
p8.quad(bottom=0, top='com_num',left='com_left', right='com_right', source=dis_hist_source, fill_color='blue', line_color='black')    

# Attach the update_plot callback to the 'value' property of select
#select1.on_change('value', update_plot)
select2.on_change('value', update_plot)
range_slider.on_change('value', update_plot)
range_slider2.on_change('value', update_plot)
range_slider3.on_change('value', update_plot)
slider.on_change('value', update_hist)

# Create layout and add to current document
layout = column(range_slider,range_slider2,range_slider3, select2, plots,slider,row(p7,p8))
#show(layout)
curdoc().clear()
curdoc().add_root(layout)
show(p1)
#cd "path_to_the_script"
#bokeh serve --show bokeh_color.py



