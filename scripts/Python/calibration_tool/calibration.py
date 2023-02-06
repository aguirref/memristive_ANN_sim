#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Created on Sun Oct 25 12:55:50 2020

@author: ngomez
"""


from dualcrossbar import DualPartitionedCrossbar, DualCrossbar
from partitionedcrossbar import PartitionedCrossbar;
from crossbarmodel import Crossbar;
from multilayercrossbar import MultiLayerCrossbar;
import scipy.io as sio;
import numpy as np;
import time;
import os.path;
import sys;

def measure_mn(syn_w):
    mn_vec=[];
    mn_vec=np.zeros(2);
    mn_vec[0]=len(syn_w);
    mn_vec[1]=len(syn_w[0]);
    return mn_vec;

try:
    directory=sys.argv[1];
except IndexError:
    directory="";


if(directory==""):
    not_calibrated=sio.loadmat('compensation_variables.mat');
else:
    completeName = os.path.join(directory, 'compensation_variables.mat') 
    not_calibrated=sio.loadmat(completeName);

matarray=not_calibrated['G_trasp'][0,0];

dimensions=matarray.ndim;

if (dimensions>2):
    number_of_layers=len(matarray);
else:
    number_of_layers=1;

syn_weights=[];
mxn=np.zeros((number_of_layers,2));

if (number_of_layers>1):
    for i in range(number_of_layers):
        syn_weights.append(np.traspose[matarray[i,0]]);
        mxn[i]=measure_mn(syn_weights[i]);
else:
    syn_weights=matarray;
    mxn=measure_mn(syn_weights);


limits=[]
limits=np.zeros(2);

images=not_calibrated['images_t10k'];
labels=not_calibrated['labels_t10k'];

limits[0]=not_calibrated['RLRS'];
limits[1]=not_calibrated['RHRS'];
line_resistance=not_calibrated['Rs'];
v_read=not_calibrated['V_read'];


if 'accuracy_calc' in not_calibrated:
    acc_calc=not_calibrated['accuracy_calc'];
else:
    acc_calc="NO";

if 'type_in' in not_calibrated:
    type_in=not_calibrated['type_in'];
else:
    type_in="SISO";
    
if 'output_side' in not_calibrated:
    output_side=not_calibrated['output_side'];
else:
    output_side="up";

if 'cal_criterion' in not_calibrated:
    cal_criterion=not_calibrated['cal_criterion'];
else:
    cal_criterion=0.05;
    
number_of_partitions=np.zeros(number_of_layers);
for i in range(number_of_layers):
    number_of_partitions[i]=not_calibrated['partitions'][i][0]*not_calibrated['partitions'][i][1];
    
if (number_of_partitions[i]>1):
    set_partitioned=1;
else:
    set_partitioned=0;

partition_m=np.zeros(number_of_layers);
partition_n=np.zeros(number_of_layers);

if(number_of_layers>1):
    for i in range(number_of_layers):
        partition_m[i]=mxn[i,0]/not_calibrated['partitions'][i,0];
        partition_n[i]=mxn[i,1]/not_calibrated['partitions'][i,1];
else:
    partition_m[0]=mxn[0]/not_calibrated['partitions'][i,0];
    partition_n[0]=mxn[1]/not_calibrated['partitions'][i,1];

if(number_of_layers>1):
    Vcal=np.zeros((2*int(mxn[0,0]+mxn[0,1])));
    Vapp=np.zeros((2*int(mxn[0,0]+mxn[0,1])));
else:
    Vcal=np.zeros((2*int(mxn[0]+mxn[1])));
    Vapp=np.zeros((2*int(mxn[0]+mxn[1])));  

for i in range(len(images)):
    for j in range(len(images[0])):
        Vcal[i]=Vcal[i]+images[i,j];
        if (type_in=='DISO'):
           Vcal[i+len(images)]=Vcal[i+len(images)]+images[i,j];     
    Vcal[i]=v_read*Vcal[i]/(len(images[0])+1);
    
crossbar_array=[];

if(set_partitioned==1):
    for i in range(number_of_layers):
        crossbar_array.append(DualPartitionedCrossbar(syn_weights,limits,"LINEAR",line_resistance,type_in,Vapp,int(number_of_partitions[i]),int(partition_m[i]),int(partition_n[i]),output_side));
    myMultiLayerCrossbar=MultiLayerCrossbar(crossbar_array,"DUAL_PARTITIONED");     
else:
    for i in range(number_of_layers):
        crossbar_array.append(DualCrossbar(syn_weights,limits,"LINEAR",line_resistance,type_in,Vapp,output_side));
    myMultiLayerCrossbar=MultiLayerCrossbar(crossbar_array,"DUAL");    

myMultiLayerCrossbar.multi_layer_calibration(Vcal,cal_criterion,type_in);

if(set_partitioned==1):
    cal_rij=[];
    weights=[];
    for i in range(number_of_layers):
        cal_rij.append(myMultiLayerCrossbar.layer_array[i].part2crossbar());
        weights.append(myMultiLayerCrossbar.layer_array[i].linear_g2w());

    dic={"G_cal":np.array(weights,dtype=np.double)}
    
    if(directory==""):
        sio.savemat('G_cal.mat',dic);
        np.savetxt("G_pos.csv", weights[0][0], delimiter=",");
        np.savetxt("G_neg.csv", weights[0][1], delimiter=",");  
    else:
        completeName = os.path.join(directory, 'G_cal.mat'); 
        sio.savemat(completeName,dic);
        completeName = os.path.join(directory, 'G_pos.csv'); 
        np.savetxt(completeName, weights[0][0], delimiter=",");
        completeName = os.path.join(directory, 'G_neg.csv');
        np.savetxt(completeName, weights[0][1], delimiter=",");

else:
    cal_rij=[];
    weights=[];
    for i in range(number_of_layers):
        cal_rij.append(myMultiLayerCrossbar.layer_array[i].cross_array[0].Rij);
        cal_rij.append(myMultiLayerCrossbar.layer_array[i].cross_array[1].Rij);
        weights.append(myMultiLayerCrossbar.layer_array[i].linear_g2w());
        
    dic={"G_cal":np.array(weights,dtype=np.double)}
    if(directory==""):
        sio.savemat('G_cal.mat',dic);
        np.savetxt("G_pos.csv", weights[0][0], delimiter=",");
        np.savetxt("G_neg.csv", weights[0][1], delimiter=",");  
    else:
        completeName = os.path.join(directory, 'G_cal.mat'); 
        sio.savemat(completeName,dic);
        completeName = os.path.join(directory, 'G_pos.csv'); 
        np.savetxt(completeName, weights[0][0], delimiter=",");
        completeName = os.path.join(directory, 'G_neg.csv');
        np.savetxt(completeName, weights[0][1], delimiter=",");

if(acc_calc=="NO"):
    sys.exit();
    
elif(acc_calc=="YES"):
    checks=0;
    imax=0;
    number_output=0;
    iout=[];
    if 'cant_muestras' in not_calibrated:
        cant_muestras=not_calibrated['cant_muestras'];
    else:
        cant_muestras=10000;
    calc_counter=0;

    for i in range(cant_muestras):
        calc_counter+=1;
        if(calc_counter%100==0):
            progress=calc_counter*100/10000;
            print("Progress: ",progress,"%\n");
            
        for j in range(len(images)):
            Vapp[j]=images[j,i];
            if(type_in=="DISO"):
                Vapp[j+len(images)]=images[j,i];
        myMultiLayerCrossbar.solve_multi_layer(Vapp);
    
        iout=myMultiLayerCrossbar.layer_array[myMultiLayerCrossbar.layers-1].get_output()
        
        imax=0;
        
        for x in range(10):
            if (iout[x]>imax):
                imax=iout[x];
                number_output=x;
        
        if (number_output==labels[0,i]):
            checks+=1;
            
    accuracy=checks/cant_muestras;
    print("Accuracy: ",accuracy*100,"%\n");
    
    
