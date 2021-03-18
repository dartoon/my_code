import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import glob
import re
import random

folders = ['233718.07+002550.6', '230402.77-003855.4', '222929.45+010438.4', '221101.45+001449.0', 
          '220642.82+003016.2', '162501.98+430931.6', '153008.91+425634.8', '150216.66+025719.8', 
          '135944.20+012809.8', '134257.16-013912.9', '133222.62+034739.9', '132441.58-015401.8', 
          '124618.51-001750.2', #'113803.73+031457.8', 
          '105458.01+043310.6', '104644.31+000329.7', 
          '101138.50+012344.6', '090654.53+021315.2', '022906.04-051428.9', '022404.85+014941.9',
          '020318.87-062321.3', '020231.14-042246.8', '020053.43-042721.8', '015141.69-000646.8', 
          '014235.36+012334.0', '013834.18-000509.3', '013736.57+005742.3', '013526.15+004135.8',
          '011227.87-003151.6', '001401.62-002600.5']
folders.sort()

ID_list, g_r_list, r_i_list = [], [], []

for ID in folders:
    folder = '../proof2close_HSC_images_5band/z_over1/' + ID + '/fit_result/'
    files = glob.glob(folder + 'fit_result_*txt')
    band = ['G', 'R', 'I', 'Z', 'Y']
    mag_list = []
    pos_list = []
    for i in range(5):
        l = [ j for j in range(len(files)) if band[i]+'-band' in files[j]]
        if l != []:
            file = files[l[0]]
            f = open(file,"r")
            string = f.read()
            lines = string.split('\n')   # Split in to \n
            l0 = [ j for j in range(len(lines)) if 'AGN mag:' in lines[j]]
            AGN_mag = lines[l0[2]]
            AGN_mag = AGN_mag.split(' ')[2:4]
            AGN_mag = [float(AGN_mag[0]), float(AGN_mag[1])]
            mag_list.append(AGN_mag)
            l1 = [ j for j in range(len(lines)) if 'PS PS' in lines[j]]
            offset = lines[l1[1]].split(' ')[-1]
            l2 = [ j for j in range(len(lines)) if 'AGN0 position:' in lines[j]]
            pos_ = lines[l2[1]]
            x0 = float((pos_.split('x: ')[1][:10]).split(' ')[0])
            x1 = float((pos_.split('x: ')[2][:10]).split(' ')[0])
            y0 = float((pos_.split('y: ')[1][:10]).split(' ')[0])
            y1 = float((pos_.split('y: ')[2][:10]).split(' ')[0])        
            pos_list.append([[x0, y0], [x1, y1]])
        else:
            mag_list.append([-99, -99])
            pos_list.append([[-99, -99],[-99, -99]])
    mag_list = np.array(mag_list)
    pos_list = np.array(pos_list)
    order_list = [[0,1]] * 5
    for i in [0, 1, 3, 4]:
        if pos_list[i][0][0] != -99: 
            dis_0 = np.sum((pos_list[i] - pos_list[2][0])**2, axis=1)
            dis_1 = np.sum((pos_list[i] - pos_list[2][1])**2, axis=1)
            if dis_0[0] > dis_0[1] and dis_1[1] > dis_1[0]:
                order_list[i] = [1, 0]
            elif dis_0[0] > dis_0[1]:
                print("Position of AGN0 werid for", ID, band[i])
            elif dis_0[0] > dis_0[1]:
                print("Position of AGN1 werid for", ID, band[i])    
    g_r_0 = mag_list[0][order_list[0][0]] - mag_list[1][0]
    r_i_0 = mag_list[1][order_list[1][0]] - mag_list[2][0]
    g_r_1 = mag_list[0][order_list[0][1]] - mag_list[1][1]
    r_i_1 = mag_list[1][order_list[1][1]] - mag_list[2][1]

    g_r_list.append([g_r_0,g_r_1])
    r_i_list.append([r_i_0, r_i_1])
    ID_list.append(ID)
#%%
from bokeh.models import ColumnDataSource, Select, RangeSlider, Slider

from bokeh.plotting import figure, output_file, show
output_file("line.html")
blue_line = np.loadtxt('Shenli_materials/stellar_blueline.txt')

p = figure(plot_width=400, plot_height=400)
p.line(blue_line[:,0], blue_line[:,1], line_width=2)
for i in range(len(g_r_list)):
    color = "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])

    qso_com = ColumnDataSource(data={
        'g_r' : g_r_list[i],
        'r_i' : r_i_list[i],
        # 'ID' : ID_list[i]
        # 'imgs' : ['/Users/Dartoon/Astro/Projects/my_code/projects/2021_dual_AGN/proof2close_HSC_images_5band/z_over1/000050.56-013055.2/fit_result/fitting2_used_aper.pdf']*2,
        'imgs' : ['/Users/Dartoon/Downloads/SED_infernece.png']*2,
    })    
    p.circle('g_r', 'r_i', source=qso_com, size=7, color=color, alpha=0.5)
from bokeh.models import HoverTool,Range1d
p.x_range = Range1d(-0.4, 1.6)
p.y_range = Range1d(-0.5, 3)
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
        </div>
        """
    )
p.add_tools(hover)
show(p)
