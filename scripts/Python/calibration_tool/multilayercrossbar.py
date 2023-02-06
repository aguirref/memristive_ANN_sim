# -*- coding: utf-8 -*-

"""
Clase que permite crear redes de múltiples capas a partir de la combinación
de crossbars

Autor: Nicolas M. Gomez

"""

"""
Bloque de importación de librerías
"""
from crossbarmodel import Crossbar
from dualcrossbar import DualCrossbar, DualPartitionedCrossbar
from partitionedcrossbar import PartitionedCrossbar
import numpy as np
import scipy.io as sio;
import math;



"""
Fin del bloque de importación de librerías
"""

"""
Bloque de declaración de clases
"""
class MultiLayerCrossbar():

#Declaración de atributos de clase

    layers=0;
    crosstype="";
    layer_array=[];
    layers_size=[];
    layers_partnum=[];
    layers_partm=[];
    layers_partn=[];
    in_size=0;

#Método de instanciación:
#Recibe:
#   *
    
#Ejecuta:
#   *
    def __init__(self,crossarray,crosstype):
        self.crosstype=crosstype;
        self.layer_array=crossarray;
        self.make_multi_layer(crossarray);

#Método que dimensiona las capas del Crossbar Multilayer
        
    def make_multi_layer(self,array):
        self.layers=len(array);
        self.layers_size=np.zeros((self.layers,2),dtype=int);
        self.layers_partnum=np.zeros((self.layers));
        self.layers_partm=np.zeros((self.layers));
        self.layers_partn=np.zeros((self.layers));
        for i in range(self.layers):
            if(self.crosstype=="DUAL"):
                self.layers_size[i,0]=array[i].m;
                self.layers_size[i,1]=array[i].n;
            if(self.crosstype=="DUAL_PARTITIONED"):
                self.layers_partnum[i]=array[i].part_number;
                self.layers_partm[i]=array[i].partm;
                self.layers_partn[i]=array[i].partn;
                self.layers_size[i,0]=array[i].m;
                self.layers_size[i,1]=array[i].n;
                     
#Método de calibración del multilayer

    def multi_layer_calibration(self,Vcal,criterion,redge):
        Vaux=[]
        iout=[]
        it_num=[];
        it_num=np.zeros((self.layers));
        for i in range(self.layers):
            print("Layer N°: ",i,"\n");
            if(self.crosstype=="DUAL"):        
                if (i==0):
                    it_num[i]=self.layer_array[i].dual_calibration(Vcal,criterion);
                    if (it_num[i]==-1):
                        return [-1, -1];
                else:
                    Vaux=[];
                    iout=[];
                    limits=[];
                    limits=np.zeros(2);
                    limits[0]=self.layer_array[i-1].LRS;
                    limits[1]=self.layer_array[i-1].HRS;
                    auxCross=DualCrossbar(self.layer_array[i-1].W,limits,"LINEAR",100e-3,redge,self.layer_array[i-1].Vapp,self.layer_array[i-1].outputside);
                    iout=auxCross.get_output();
                    Vaux=np.zeros(2*(self.layers_size[i,0]+self.layers_size[i,1]))
                    for j in range (self.layers_size[i-1,1]):
                        Vaux[j]=0.5*self.sigmoid(1e6*iout[j]);
                    it_num[i]=self.layer_array[i].dual_calibration(Vaux,criterion);
                    del auxCross;
            elif(self.crosstype=="DUAL_PARTITIONED"):
                if(i==0):
                    it_num[i]=self.layer_array[i].partitioned_dual_calibration(Vcal,criterion);
                    if (it_num[i]==-1):
                        return[-1,-1];
                else:
                    Vaux=[];
                    iout=[];
                    limits=[];
                    limits=np.zeros(2);
                    limits[0]=self.layer_array[i-1].LRS;
                    limits[1]=self.layer_array[i-1].HRS;
                    auxCross=DualPartitionedCrossbar(self.layer_array[i-1].W,limits,"LINEAR",100e-3,redge,self.layer_array[i-1].Vapp,self.layer_array[i-1].part_number,self.layer_array[i-1].partm,self.layer_array[i-1].partn,self.layer_array[i-1].outputside);
                    iout=auxCross.get_output();
                    Vaux=np.zeros(2*(self.layers_size[i,0]+self.layers_size[i,1]))
                    for j in range (self.layers_size[i-1,1]):
                        Vaux[j]=0.5*self.sigmoid(1e6*iout[j]);
                    it_num[i]=self.layer_array[i].partitioned_dual_calibration(Vaux,criterion);
                    del auxCross;
        return it_num;             
    
#Método de resolución

    def solve_multi_layer(self,vapp):
        Vaux=[];
        iout=[];
        for i in range(self.layers):
            if(self.crosstype=="DUAL"):
                if (i==0):
                    self.layer_array[i].solve_dual(vapp);
                else:
                    iout=[];
                    Vaux=[];
                    iout=np.zeros(self.layers_size[i-1,1]);
                    Vaux=np.zeros(2*(self.layers_size[i,0]+self.layers_size[i,1]))
                    iout=self.layer_array[i-1].get_output();
                    for j in range (self.layers_size[i-1,1]):
                        Vaux[j]=0.5*self.sigmoid(1e6*iout[j]);
                    self.layer_array[i].solve_dual(Vaux);
            if(self.crosstype=="DUAL_PARTITIONED"):
                if (i==0):
                    self.layer_array[i].solve_partitioned_dual(vapp);
                else:
                    iout=[];
                    Vaux=[];
                    iout=np.zeros(self.layers_size[i-1,1]);
                    Vaux=np.zeros(2*(self.layers_size[i,0]+self.layers_size[i,1]))
                    iout=self.layer_array[i-1].get_output();
                    for j in range (self.layers_size[i-1,1]):
                        Vaux[j]=0.5*self.sigmoid(1e6*iout[j]);
                    self.layer_array[i].solve_partitioned_dual(Vaux);     

#Función sigmoidal    

    def sigmoid(self,x):
        return 1 / (1 + math.exp(-x))
        
"""
Fin del bloque de declaración de clases
"""

"""
Bloque de declaración de métodos
"""

    
"""
Fin del bloque de declaración de métodos
"""


"""
Bloque principal de código
"""

####################TEST###########################
#TEST ROUTINE

# imageset=sio.loadmat('Test/multilayer/3_layers/64_54_10/Scaled_Conjugate_Gradient/train_1/G_values.mat');

# matarray=imageset['G_real'];
# imagenes=imageset['images_test'];
# labels=imageset['labels_test']

# syn_weights1=np.transpose(matarray[0,0]);
# syn_weights2=np.transpose(matarray[1,0]);

# limites=[100e3,1e6];

# Vcal=np.zeros((2*(64+54)));
# Vapp1=np.zeros((2*(64+54)));
# Vapp2=np.zeros((2*(54+10)));



# myFirstCrossbar = DualCrossbar(syn_weights1,limites,"LINEAR",20,"SISO",Vapp1);
# mySecondCrossbar = DualCrossbar(syn_weights2,limites,"LINEAR",20,"SISO",Vapp2);

# array=[];
# array.append(myFirstCrossbar);
# array.append(mySecondCrossbar)

# myMultiLayer = MultiLayerCrossbar(array, "DUAL");

# checks=0;
# imax=0;
# number_output=0;
# iout=[];

# for i in range(len(imagenes)):
#     for j in range(len(imagenes[0])):
#         Vcal[i]=Vcal[i]+imagenes[i,j];
#     Vcal[i]=Vcal[i]/(len(imagenes[0])+1);

# itnum=[];
# itnum=myMultiLayer.multi_layer_calibration(Vcal);

# for i in range (100):
#     Vapp=np.zeros((2*(64+54)));
#     for j in range (64):
#         Vapp[j]=imagenes[j,i];
#         #Vapp[j+64]=imagenes[j,i];
#     myMultiLayer.solve_multi_layer(Vapp);

#     iout=myMultiLayer.layer_array[myMultiLayer.layers-1].get_output();
    
#     imax=0;
#     for x in range (10):
#         if (iout[x]>imax):
#             imax=iout[x];
#             number_output=x;
    
#     if (number_output==labels[0,i]):
#         checks+=1;

# print(checks);
        

        
    
