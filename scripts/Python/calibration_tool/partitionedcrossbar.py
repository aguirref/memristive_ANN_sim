# -*- coding: utf-8 -*-
"""
Clase que permite crear un Crossbar particionado en p particiones de mxn.
Para la resolución de las particiones, se utiliza la clase Crossbar.

Autor: Nicolas M. Gomez

"""

"""
Bloque de importación de librerías
"""
from crossbarmodel import Crossbar
import numpy as np

"""
Fin del bloque de importación de librerías
"""

"""
Bloque de declaración de clases
"""
class PartitionedCrossbar():

#Declaración de atributos de clase
    
    partitions=[]
    part_Vapp=[]
    partnum=0;
    partperrow=0;
    partpercolumn=0;
    partm=0;
    partn=0;
    m=0;
    n=0;
    Rij=[]
    Vapp=[]
    cross_array=[]
    outputside="up";

#Método de instanciación:
#Recibe:
#   *Matriz Rij
#   *Valor de las R de linea "internas"
#   *String que identifica al tipo de entrada/salida
#   *Vector de tensiones aplicadas
#   *Número de particiones
#   *Filas por partición
#   *Columnas por partición
    
#Ejecuta:
#   *Método para estructurar las particiones y el array correspondiente
#   *Método para la construcción de las particiones Rij
#   *Método para la resolución de las p particiones
    
    def __init__(self,vecrij,limits,valorrl,redge,vapp,partnumber,partm,partn,outside):
        self.partnum=partnumber;
        self.Rij=vecrij;
        self.Vapp=vapp;
        self.outputside=outside;
        self.set_partitions(vecrij,partnumber,partm,partn);
        self.make_partitions();
        self.solve_partitions(limits,valorrl,redge);
        
#Método que estructura las particiones y el vector que las agrupa
        
    def set_partitions(self,Rij,pnumber,rows,col):
        self.m=len(Rij);
        self.n=len(Rij[0]);
        self.partm=rows;
        self.partn=col;

        self.partitions=np.zeros((pnumber,rows,col));
        self.part_Vapp=np.zeros((pnumber,2*self.partm+2*self.partn,1));
        
        self.partperrow=int(self.n/col);
        self.partpercolumn=int(self.m/rows)

#Construcción de las p particiones de Rij y Vapp
        
    def make_partitions(self):

        
        partcounter=0;
        
        for j in range(self.partpercolumn):
            for i in range(self.partperrow):
                for x in range(self.partm):
                    for y in range(self.partn):
                        self.partitions[partcounter,x,y]=self.Rij[j*self.partm+x,i*self.partn+y];
                
                for z in range(2*(self.partm+self.partn)):
                    if (z<self.partm):
                        self.part_Vapp[partcounter,z,0]=self.Vapp[j*self.partm+z];
                    elif (z>=self.partm and z<2*self.partm):
                        self.part_Vapp[partcounter,z,0]=self.Vapp[self.m+j*self.partm+z-self.partm];
                    elif (z>=2*self.partm and z<(2*self.partm+self.partn)):
                        self.part_Vapp[partcounter,z,0]=self.Vapp[2*self.m+i*self.partn+z-2*self.partm];
                    else:
                        self.part_Vapp[partcounter,z,0]=self.Vapp[2*self.m+self.n+i*self.partn+z-2*self.partm-self.partn];
                partcounter+=1;
   
#Método que crea y resuelve los p Crossbar
                
    def solve_partitions(self,limits,valorrl,redge):
        self.flush_partitioned();
        for i in range(self.partnum):            
            self.cross_array.append(Crossbar(self.partitions[i,:,:],limits,valorrl,redge,self.part_Vapp[i,:,:],self.outputside));
            
#Método utilizado para resolver con nuevas vapp
                
    def solve_partitioned(self,vapp):
        self.Vapp=vapp;
        partcounter=0;
        
        for j in range(self.partpercolumn):
            for i in range(self.partperrow):
             
                for z in range(2*(self.partm+self.partn)):
                    if (z<self.partm):
                        self.part_Vapp[partcounter,z,0]=self.Vapp[j*self.partm+z];
                    elif (z>=self.partm and z<2*self.partm):
                        self.part_Vapp[partcounter,z,0]=self.Vapp[self.m+j*self.partm+z-self.partm];
                    elif (z>=2*self.partm and z<(2*self.partm+self.partn)):
                        self.part_Vapp[partcounter,z,0]=self.Vapp[2*self.m+i*self.partn+z-2*self.partm];
                    else:
                        self.part_Vapp[partcounter,z,0]=self.Vapp[2*self.m+self.n+i*self.partn+z-2*self.partm-self.partn];
                partcounter+=1;
        
        for i in range(self.partnum):
            self.cross_array[i].cross_solve(self.part_Vapp[i,:,:]);
            
#Método de calibración del crossbar particionado
    def partitioned_calibration(self,Vcal,criterion):
        it_num=[];
        it_num=np.zeros((self.partnum));
        itmax=0;
        
        self.Vapp=Vcal;
        partcounter=0;
        
        for j in range(self.partpercolumn):
            for i in range(self.partperrow):
             
                for z in range(2*(self.partm+self.partn)):
                    if (z<self.partm):
                        self.part_Vapp[partcounter,z,0]=self.Vapp[j*self.partm+z];
                    elif (z>=self.partm and z<2*self.partm):
                        self.part_Vapp[partcounter,z,0]=self.Vapp[self.m+j*self.partm+z-self.partm];
                    elif (z>=2*self.partm and z<(2*self.partm+self.partn)):
                        self.part_Vapp[partcounter,z,0]=self.Vapp[2*self.m+i*self.partn+z-2*self.partm];
                    else:
                        self.part_Vapp[partcounter,z,0]=self.Vapp[2*self.m+self.n+i*self.partn+z-2*self.partm-self.partn];
                partcounter+=1;

        for i in range(self.partnum):
            it_num[i]=self.cross_array[i].calibrate(self.part_Vapp[i,:,:],criterion);
            if(it_num[i]==-1):
                itmax=-1;
            elif(it_num[i]>itmax):
                itmax=it_num[i];
        return itmax;

#Método utilizado para obtener las corrientes de salida
            
    def get_output(self):
        out=[];
        out=np.zeros(self.n);
        
        for i in range (self.partperrow):
            aux=np.zeros(self.partn);
            for j in range (self.partpercolumn):
                aux=aux+self.cross_array[i+j*self.partperrow].get_output();
            for x in range(self.partn):       
                out[x+i*self.partn]=aux[x];
                
        return out;
    
    def flush_partitioned(self):
        self.cross_array=[];
    
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

# f=16;
# c=16;

# Rij=[]
# Rij=np.zeros((f,c))

# for i in range (f):
#     for j in range (c):
#         if (i==j):
#             Rij[i,j]=100e3;
#         else:
#             Rij[i,j]=10e6;
            
# Vapp=np.zeros((2*(f+c)));
# Vapp[3]=0.5;
# # Vapp[19]=0.5;
# # Vapp[15]=0.5;
# # Vapp[7]=0.5;
# #Vapp[32]=0.5;
# #Vapp[48]=0.5;

            
# mypartcrossbar = PartitionedCrossbar(Rij,100e-3,"SISO",Vapp,16,4,4);

# io=[];

# io=mypartcrossbar.get_output();