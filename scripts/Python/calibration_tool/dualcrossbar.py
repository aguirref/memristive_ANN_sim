# -*- coding: utf-8 -*-
"""
Clase que permite crear un doble Crossbar en función de los pesos sinápticos
recibidos como argumento. En un Crossbar se resuelve el modelo de memristores
que identifican a los pesos positivos. En el otro Crossbar, se realiza el
mismo procedimiento para los negativos.


Autor: Nicolas M. Gomez

"""

"""
Bloque de importación de librerías
"""
from crossbarmodel import Crossbar
from partitionedcrossbar import PartitionedCrossbar
import numpy as np
import multiprocessing as mp


"""
Fin del bloque de importación de librerías
"""

"""
Bloque de declaración de clases
"""
class DualCrossbar():

#Declaración de atributos de clase

    m=0;
    n=0;
    Rij=[];
    Vapp=[];
    cross_array=[];
    W=[];
    HRS=0;
    LRS=0;
    VlineWL=[];
    VlineBL=[];
    outputside="up";

#Método de instanciación:
#Recibe:
#   *Matriz de pesos
#   *Vector que indica HRS y LRS (en ese orden)
#   *String que identifica el modo de conversión de pesos sinápticos a
#    resistencias (actualmente, solo linear disponible)
#   *Valor de resistencias de línea internas
#   *String que identifica al tipo de entrada/salida
#   *Vector de tensiones aplicadas
    
#Ejecuta:
#   *Método para dimensionar el problema a resolver
#   *Método para la construcción de las dos matrices Rij
#   *Método para la resolución de los Crossbars
    
    def __init__(self,w,limits,mode,valorrl,redge,vapp,outside):

        self.W=w;
        self.Vapp=vapp;
        self.LRS=limits[0];
        self.HRS=limits[1];
        self.outputside=outside;
        
        self.get_row_column();
        if (mode=='LINEAR'):
            self.linear_w2g();
        self.create_crossbar(limits,valorrl,redge);
        
#Método que dimensiona la matriz de pesos
            
    def get_row_column(self):
        self.m=len(self.W);
        self.n=len(self.W[0]);
        
        self.Rij=np.zeros((2,self.m,self.n));
        
        
#Método crea las matrices Rij a partir de la matriz de pesos
                
    def linear_w2g(self):
        
        max = self.W[0,0];
        aux = 0;
        
        Gmin=1/self.HRS;
        Gmax=1/self.LRS;
        
        for i in range (self.m):
            for j in range (self.n):
                if (max<abs(self.W[i,j])):
                    max=abs(self.W[i,j]);
        
        for i in range (self.m):
            for j in range (self.n):
                aux=self.W[i,j]/max;
                if (aux>0):
                    self.Rij[0,i,j]=1/(aux*(Gmax-Gmin)+Gmin);
                    self.Rij[1,i,j]=1/(Gmin);
                else:
                    self.Rij[0,i,j]=1/(Gmin);
                    self.Rij[1,i,j]=1/(np.absolute(aux)*(Gmax-Gmin)+Gmin);
                
#Método que devuelve el tensor de pesos a partir de las matrices Rij

    def linear_g2w(self):
        
        Gmin=1/self.HRS;
        Gmax=1/self.LRS;
        
        G=np.zeros((2,self.m,self.n));
        
        for i in range (self.m):
            for j in range (self.n):
                G[0,i,j]=((1/self.Rij[0,i,j])-Gmin)/(Gmax-Gmin);
                G[1,i,j]=((1/self.Rij[1,i,j])-Gmin)/(Gmax-Gmin);

        return G;

#Crea y resuelve los elementos Crossbar del array

    def create_crossbar(self,limits,valorrl,redge):
        self.flush_dual();
        self.cross_array.insert(0,Crossbar(self.Rij[0,:,:],limits,valorrl,redge,self.Vapp,self.outputside));
        self.cross_array.insert(1,Crossbar(self.Rij[1,:,:],limits,valorrl,redge,self.Vapp,self.outputside));
        self.VlineBL=self.cross_array[0].VlineBL-self.cross_array[1].VlineBL;
        self.VlineWL=self.cross_array[0].VlineWL-self.cross_array[1].VlineWL;
        
#Solución al dual Crossbar
                    
    def solve_dual(self,vapp):
        self.Vapp=vapp;
        self.cross_array[0].cross_solve(vapp);
        self.cross_array[1].cross_solve(vapp);
        self.VlineBL=self.cross_array[0].VlineBL-self.cross_array[1].VlineBL;
        self.VlineWL=self.cross_array[0].VlineWL-self.cross_array[1].VlineWL;

#Calibración del Crossbar Dual
    
    def dual_calibration(self, Vcal,criterion):
        self.Vapp=Vcal;

        print("Calibración de pesos positivos\n");
        it_num1=self.cross_array[0].calibrate(Vcal,criterion);
        print("Calibración de pesos negativos\n");
        it_num2=self.cross_array[1].calibrate(Vcal,criterion);
        
        
        if(it_num1>=it_num2):
            return it_num1;
        else:
            return it_num2;


#Método utilizado para obtener las corrientes de salida
        
    def get_output(self):
        out=[];
        out=np.zeros(self.n);
        if(self.outputside=="down"):
            for i in range(self.n):
                out[i]=(self.VlineBL[self.m-1,i]-self.Vapp[2*self.m+self.n+i])/self.cross_array[0].Rbl2
        else:
            for i in range(self.n):
                out[i]=(self.VlineBL[0,i]-self.Vapp[2*self.m+i])/self.cross_array[0].Rbl1            
        
        return out;
    
#Método que limpia los elementos del array
        
    def flush_dual(self):
        self.cross_array=[];
    
#A continuación se declara la clase DualPartionedCrossbar, la cual hereda
#la clase DualCrossbar y, a su vez, implementa la clase PartitionedCrossbar
#para simular el dispositivo dual particionado.
        
class DualPartitionedCrossbar(DualCrossbar):
    
    partitioned_cross_array=[];
    part_number=0;
    partm=0;
    partn=0;

#Definición de método de inicialización, al cual, además de los parámetros necesarios,
#para instanciar DualCrossbar, se le deben pasar también aquellos necesarios para generar las particiones:
#   *Número de particiones
#   *Filas de partición
#   *Columnas de partición

    def __init__(self,w,limits,mode,valorrl,redge,vapp,partnumber,partm,partn,outside):
        super().__init__(w,limits,mode,valorrl,redge,vapp,outside);
        self.part_number=partnumber;
        self.partm=partm;
        self.partn=partn;
        self.create_partitioned_dual(limits,valorrl,redge,partnumber,partm,partn,outside);

#Método utilizado para la resolución del Crossbar Dual Particionado
        
    def create_partitioned_dual(self,limits,valorrl,redge,pnumber,pm,pn,outside):
        self.flush_dual();
        self.partitioned_cross_array.append(PartitionedCrossbar(self.Rij[0,:,:],limits,valorrl,redge,self.Vapp,pnumber,pm,pn,outside))
        self.partitioned_cross_array.append(PartitionedCrossbar(self.Rij[1,:,:],limits,valorrl,redge,self.Vapp,pnumber,pm,pn,outside))

#Método utilizado para la resolución del Crossbar Dual Particionado
        
    def solve_partitioned_dual(self,vapp):
        self.partitioned_cross_array[0].solve_partitioned(vapp);
        self.partitioned_cross_array[1].solve_partitioned(vapp);

#Método de calibración del crossbar dual particionado
    
    def partitioned_dual_calibration(self,Vcal,criterion):
        self.Vapp=Vcal;

        print("Calibración de pesos positivos\n");
        it_num1=self.partitioned_cross_array[0].partitioned_calibration(Vcal,criterion);
        print("Calibración de pesos negativos\n");
        it_num2=self.partitioned_cross_array[1].partitioned_calibration(Vcal,criterion);
        
        if(it_num1>=it_num2):
            return it_num1;
        else:
            return it_num2;
        
#Método que rearma los crossbar

    def part2crossbar(self):
        aux=[];
        aux=np.zeros((2,self.m,self.n));
        partsperrow=int(self.n/self.partn);
        partspercolumn=int(self.m/self.partm);
        partcounter=0;
        cross_row=0;
        cross_column=0;
       
        for i in range(partspercolumn):
            for j in range(partsperrow):
                for x in range(self.partn):
                    for y in range(self.partm):
                        aux[0,cross_row+y,cross_column+x]=self.partitioned_cross_array[0].cross_array[partcounter].Rij[y,x];
                cross_column=cross_column+self.partn;
                partcounter+=1;
            cross_row=cross_row+self.partm;
            cross_column=0;
            
        partcounter=0;
        cross_row=0;
        cross_column=0;
        
        for i in range(partspercolumn):
            for j in range(partsperrow):
                for x in range(self.partn):
                    for y in range(self.partm):
                        aux[1,cross_row+y,cross_column+x]=self.partitioned_cross_array[1].cross_array[partcounter].Rij[y,x];
                cross_column=cross_column+self.partn;
                partcounter+=1;
            cross_row=cross_row+self.partm;
            cross_column=0;
            
        self.Rij=aux;
            
        return aux;
        
                                  
            
                

#Método que presenta las corrientes de salida
       
    def get_output(self):
        out=[];
        out=np.zeros(self.n);
        
        out=(self.partitioned_cross_array[0].get_output())-(self.partitioned_cross_array[1].get_output());
        
        return out;
    
    def flush_dual(self):
        self.partitioned_cross_array=[];
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

# f=64;
# c=10;

# Rij=[]
# Rij=np.zeros((f,c))

# for i in range (f):
#     for j in range (c):
#         if (i==j):
#             Rij[i,j]=1;
#         else:
#             Rij[i,j]=-1;
            
# limits = [100e3,10e6];

# Vapp=np.zeros((2*(f+c)));

# #Vapp[14]=0.5;
# Vapp[3]=0.5;
# Vapp[5]=0.5;

# myDual = DualCrossbar(Rij, limits, "LINEAR", 100e-3, "SISO", Vapp);



# Vapp=np.zeros((2*(74)));

# for j in range (64):
#     Vapp[j]=imagenes[j,2];
    
# myCrossbar = DualPartitionedCrossbar(rij,limites,"LINEAR",100e-3,"SISO",Vapp,8,16,5);

# iout=[]
# iout=myCrossbar.get_output();
            
# mydual = DualPartitionedCrossbar(Rij,limits,"LINEAR",100e-3,"SISO",Vapp,16,4,4);

# iout=[]

# iout=mydual.get_output();