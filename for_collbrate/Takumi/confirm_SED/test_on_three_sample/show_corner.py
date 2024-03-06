#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 09:54:48 2023

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

#%%
import pickle, glob
file = glob.glob('./esti_smass_lines/112911/'+'chain_*_phys.cpkl')
value = pickle.load(open(file[0],'rb'))
import seaborn
import pandas as pd 

prop = list(value['chain'].keys())
for i in range(len(prop)):
    if i == 1:
        value['chain'][prop[i]] = 10**value['chain'][prop[i]]

df = pd.DataFrame(data=value['chain'])

g = seaborn.pairplot(df,corner=True, diag_kind="kde", plot_kws=dict(marker="+", linewidth=1))
# g.map_lower(seaborn.kdeplot, levels=4, color=".2")

values = []
for key in df.keys():
    values.append([np.percentile(df[key], 16), np.median(df[key]), np.percentile(df[key], 84)])
for i in range(len(g.axes.ravel())):
    ax = g.axes.ravel()[i]
    if ax != None:
        ax.spines['top'].set_visible(True) 
        ax.spines['right'].set_visible(True)
        ax.spines['left'].set_visible(True)
        j = i % len(values)
        ax.axvline(x=values[j][1], ls='--', linewidth=1.6, c='coral')
        ax.axvline(x=values[j][0], ls='--', linewidth=1.6, c='coral')
        ax.axvline(x=values[j][2], ls='--', linewidth=1.6, c='coral')
        if i %5 == 0:
            title_fmt=".2f"
            fmt = "{{0:{0}}}".format(title_fmt).format
            title = r"${{{0}}}_{{-{1}}}^{{+{2}}}$"
            title = title.format(fmt(values[j][1]), fmt(values[j][1]-values[j][0]), fmt(values[j][2]-values[j][1]))
            ax.set_title(title, fontsize=16)
            # ax.set_title(df.keys()[i/5]+r'$\displaystyle\substack{1\2}$'.format(values[j][1], values[j][2]-values[j][1], values[j][1]-values[j][0]))
            # ax.set_title(r'\frac{-e^{i\pi}}{2^n}$!')
            # ax.set_title(r'\frac{-e^{i\pi}}{2^n}$!', fontsize=16, color='r')
            # ax.set_title(r'\TeX\ is Number $\displaystyle\sum_{n=1}^\infty'
            #              r'\frac{-e^{i\pi}}{2^n}$!', fontsize=16, color='r')
            # ax.set_title(df.keys()[i/5]+r'={0}$\pm$')

