#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 10:50:18 2017

@author: dxh
"""
import numpy as np
import copy
def roundme(dicts, prec=3):
    if isinstance(dicts,float):
        dicts= round(dicts, prec)
        return dicts
    if isinstance(dicts, np.ndarray):
        dicts=np.around(dicts,prec)
        return dicts
    if isinstance(dicts, dict):
        for k, v in dicts.items():
            if isinstance(v,float):
                dicts[k] = round(v, prec)
            if isinstance(v,np.ndarray):
                for i in range(len(v)):
                    v[i]=round(v[i], prec)
        return dicts
    if isinstance(dicts, tuple):
        for i in range(len(dicts)):
            for k, v in dicts[i].items():
                dicts[i][k] = round(v, prec)
        return dicts
    else:
        return dicts

def roundme_str(dicts, prec=5):
    if isinstance(dicts,float):
        dicts= '{0:.{1}f}'.format(dicts,prec) #round(dicts, prec)
        return dicts
    if isinstance(dicts, np.ndarray):
        str_list = []
        for i in range(len(dicts)): 
            str_list.append('{0:.{1}f}'.format(dicts[i],prec) ) 
        return str_list
    if isinstance(dicts, dict):
        dicts_ = copy.deepcopy(dicts)
        for k, v in dicts.items():
            if isinstance(v,float):
                dicts_[k] = '{0:.{1}f}'.format(v,prec) #round(v, prec)
            if isinstance(v, np.ndarray):
                str_list = []
                for i in range(len(v)):
                    str_list.append('{0:.{1}f}'.format(v[i],prec) ) 
                dicts_[k] = str_list
        return dicts_
    if isinstance(dicts, tuple):
        dicts_ = copy.deepcopy(dicts)
        for i in range(len(dicts_)):
            for k, v in dicts[i].items():
                dicts_[i][k] = '{0:.{1}f}'.format(v,prec)  #round(v, prec)
        return dicts
    if isinstance(dicts, list):
        if isinstance(dicts[0], float):
            str_list = []
            for i in range(len(dicts)):
                str_list.append('{0:.{1}f}'.format(dicts[i],prec) ) 
            return str_list  
        elif isinstance(dicts[0], dict):
            str_list = []
            for i in range(len(dicts)):
                dicts_ = copy.deepcopy(dicts[i])
                for k, v in dicts[i].items():
                    if isinstance(v,float):
                        dicts_[k] = '{0:.{1}f}'.format(v,prec) #round(v, prec)
                    if isinstance(v, np.ndarray):
                        str_list = []
                        for i in range(len(v)):
                            str_list.append('{0:.{1}f}'.format(v[i],prec) ) 
                        dicts_[k] = str_list   
                str_list.append(dicts_)
            return str_list  
    else:
        return dicts
