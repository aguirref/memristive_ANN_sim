#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 10:37:15 2020

@author: ngomez
"""

import scipy.io as sio;
import numpy as np;

not_calibrated=sio.loadmat("compensation_variables.mat");

matarray=not_calibrated['G'][0,0];

calibrated_dict=sio.loadmat("G_cal.mat");

dictlist=[];

for key, value in calibrated_dict.items():
    temp = [key,value]
    dictlist.append(temp)

calibrated=value


G_tot_cal=calibrated[0,0,:,:]-calibrated[0,1,:,:];


comparison=matarray-G_tot_cal;


maximo=np.max(abs(comparison));
