# -*- coding: utf-8 -*-
"""
Modelización de un crossbar que contempla la influencia de las resistencias de
línea, basada en las ecuaciones matriciales planteadas en el paper 
"A Comprehensive Crossbar Array Model With Solutions for Line Resistance and
Nonlinear Device Characteristics " de An Chen

Autor: Nicolas M. Gomez

"""

"""
Bloque de importación de librerías
"""
import numpy as np
import scipy.sparse.linalg as sp
#from repoze.lru import lru_cache
from random import randint
import time

"""
Fin del bloque de importación de librerías
"""

"""
Bloque de declaración de clases
"""
class Crossbar:
    Rline=0
    Rij=[];
    A=[];
    B=[];
    C=[];
    D=[];
    ABCD=[];
    Vline=[];
    E=[];
    ABCD_inv=[];
    Rwl1=0;
    Rwl2=0;
    Rbl1=0;
    Rbl2=0;
    VlineWL=[];
    VlineBL=[];
    hliner=10e6;
    lliner=1e-6;
    Vapp=[];
    cal=[];
    HRS=0;
    LRS=0;
    outputside="up";

#Inicializador de cada instancia de clase. Actualmente:
#           *Se pasan como parámetros: vector de elementos Rij,
#            valor de la resistencia de linea, string que permite setear
#            las resistencias de linea de los bordes y el
#            vector de tensiones aplicadas     
#           *Se miden las dimensiones de Rij para pasar el número de filas y columnas
#            atributos de clase correspondientes
#           *Se llama a los métodos para construir las matrices a, b, c y d
#            individualmente y luego se conforma la matriz de coeficientes
#           *Se conforman las matrices E y V y se ejecuta el método resolutor
#           *El resultado quedará almacenado en Vline[]
#           *Se incorpora un método para mapear las tensiones de nodo en formato matricial
    
    def __init__(self,vecrij,limits,valorrl,redge,vapp,outside):
        self.get_row_column(vecrij);
        self.outputside=outside;
        self.Rline=valorrl
        self.LRS=limits[0];
        self.HRS=limits[1];
        self.lliner=valorrl;
        self.Rij=vecrij;
        self.Vapp=vapp;
        self.set_edge_rline(redge);
        self.matrix_a();
        self.matrix_b();
        self.matrix_c();
        self.matrix_d();
        self.matrix_abcd();
        self.matrix_vline();
        self.matrix_E(vapp)
        self.matrix_solve();
        self.vline_mapping();

#Método que obtiene el número de columnas y filas a partir de la matriz Rij
        
    def get_row_column(self, Rij):
        self.f=len(Rij);
        self.c=len(Rij[0]);

    
#Método que llena a la matriz Rij con un solo elemento en su diagonal
        
    def set_rij_diag1element(self,Rij,valor):
        n=self.c;
        m=self.f;
        for i in range(m):
            for j in range(n):
                if(i==j):
                    Rij[i,j]=valor;
                else:
                    Rij[i,j]=1e6;
  
#Método que permite determinar Rij con valores aleatorios  
                
    def set_rij_random(self,Rij,a,b):
        n=self.c;
        m=self.f;
        for i in range(n):
            for j in range(m):
                Rij[i,j]=randint(a,b);


#Método utilizado para setear las resistencia de línea de los bordes
#Recibe como parámetro un stringe que responderá a cualquiera de las siguientes
#opciones:
#       *SISO: Single input, single output
#       *DISO: Double input, single output
#       *SIDO: Single input, double output
#       *DIDO: Double input, double output
                
    def set_edge_rline(self, redge):
        if(redge=="SISO"):
            if(self.outputside=="down"): 
                self.Rwl1=self.lliner;
                self.Rwl2=self.hliner;
                self.Rbl1=self.hliner;
                self.Rbl2=self.lliner;
            else:
                self.Rwl1=self.lliner;
                self.Rwl2=self.hliner;
                self.Rbl1=self.lliner;
                self.Rbl2=self.hliner;
        elif(redge=="DISO"):
            if(self.outputside=="down"): 
                self.Rwl1=self.lliner;
                self.Rwl2=self.lliner;
                self.Rbl1=self.hliner;
                self.Rbl2=self.lliner;
            else:
                self.Rwl1=self.lliner;
                self.Rwl2=self.lliner;
                self.Rbl1=self.lliner;
                self.Rbl2=self.hliner;
        elif(redge=="DIDO"):
            self.Rwl1=self.lliner;
            self.Rwl2=self.lliner;
            self.Rbl1=self.lliner;
            self.Rbl2=self.lliner;
        elif(redge=="SIDO"):
            self.Rwl1=self.lliner;
            self.Rwl2=self.hliner;
            self.Rbl1=self.lliner;
            self.Rbl2=self.lliner;
        else:
            print("\nAl instanciar la clase, el parámetro debe ser un string (ver documentación de clase)\n")

            
#Método que consume la propia clase para construir la matriz A  

    def matrix_a(self):
        n=self.c;
        m=self.f;
        self.A=np.zeros((m*n,m*n))
        
        constant=2/self.Rline;
        edgewl1constant=(1/self.Rline)+(1/self.Rwl1);
        edgewl2constant=(1/self.Rline)+(1/self.Rwl2);
        
        for i in range (m):
            if (m>=n):
                for aux in range (m):
                    for j in range (n):
                        if (aux==j):
                            if(aux==0):
                                self.A[aux+(i*n),j+(i*n)]=(1/self.Rij[i,j])+edgewl1constant;
                            elif(aux==(n-1)):
                                self.A[aux+(i*n),j+(i*n)]=(1/self.Rij[i,j])+edgewl2constant;
                            else:
                                self.A[aux+(i*n),j+(i*n)]=(1/self.Rij[i,j])+constant;
                        if(((aux+1)==j) or ((aux-1)==j)):
                            if(((j+(i*n))<(n*m)) and ((aux+(i*n))<(n*m))):
                                if(0<=aux<n):
                                    self.A[aux+(i*n),j+(i*n)]=-(1/self.Rline);
                                
            else:   
                for aux in range (n):
                    for j in range (n):
                        if (aux==j):
                            if(aux==0):
                                self.A[aux+(i*n),j+(i*n)]=(1/self.Rij[i,j])+edgewl1constant;
                            elif(aux==(n-1)):
                                self.A[aux+(i*n),j+(i*n)]=(1/self.Rij[i,j])+edgewl2constant;
                            else:
                                self.A[aux+(i*n),j+(i*n)]=(1/self.Rij[i,j])+constant;
                        if((aux+1)==j or (aux-1)==j):
                            if(((j+(i*n))<(n*m)) and ((aux+(i*n))<(n*m))):
                                if(0<=aux<n):
                                    self.A[aux+(i*n),j+(i*n)]=-(1/self.Rline);
                                
                                
#Método que consume la propia clase para construir la matriz B  
       
    def matrix_b(self):
        n=self.c;
        m=self.f;
        self.B=np.zeros((m*n,m*n))
        
        for i in range (m): 
            if (m>=n):
                for aux in range (m):
                    for j in range (n):
                        if (aux==j):
                            self.B[aux+(i*n),j+(i*n)]=-1/self.Rij[i,j];
            else:
                for aux in range (n):
                    for j in range (n):
                        if (aux==j):
                            self.B[aux+(i*n),j+(i*n)]=-1/self.Rij[i,j];
    
    


#Método que consume la propia clase para construir la matriz C  

    def matrix_c(self):
        n=self.c;
        m=self.f;
        self.C=np.zeros((m*n,m*n))

        for j in range(n):
            for i in range(m):
                if(i<m and j<n):
                    self.C[(m*j)+i,n*i+j]=1/(self.Rij[i,j])

            
#Método que consume la propia clase para construir la matriz D  
       
    def matrix_d(self):
        n=self.c;
        m=self.f;
        self.D=np.zeros((m*n,m*n));
        
        constant=2/self.Rline;
        edgebl1constant=(1/self.Rline)+(1/self.Rbl1);
        edgebl2constant=(1/self.Rline)+(1/self.Rbl2);

        for j in range(n):
            for i in range(m):
                if(i<m and j<n):
                    if(i==0):
                        self.D[(m*j)+i,n*i+j]=-((1/(self.Rij[(i,j)]))+edgebl1constant);
                    elif (i==(m-1)):
                        self.D[(m*j)+i,n*i+j]=-((1/(self.Rij[(i,j)]))+edgebl2constant);
                    else:
                        self.D[(m*j)+i,n*i+j]=-((1/(self.Rij[(i,j)]))+constant);
                    if(i<m-1):
                        self.D[(m*j)+i,n*(i+1)+j]=1/self.Rline;
                    if(i>0):
                        self.D[(m*j)+i,n*(i-1)+j]=1/self.Rline;

                           
#Método que consume la propia clase para componer ABCD
                       
    def matrix_abcd(self):
        n=self.c;
        m=self.f;
        self.ABCD=np.zeros((2*m*n,2*m*n));
        
        for i in range(2*m*n):
            if(i<2*m*n):
                if(i<m*n):
                    for j in range(m*n):
                        if (j<m*n):
                            self.ABCD[i,j]=self.A[i,j];
                            self.ABCD[i,j+(m*n)]=self.B[i,j];
                else:
                    for j in range(m*n):
                        if (j<m*n):
                            self.ABCD[i,j]=self.C[i-m*n,j];
                            self.ABCD[i,j+m*n]=self.D[i-m*n,j];

#Método llamado por la clase para declarar un vector de 2*m*n lleno de ceros
#en el que se almacenarán los valores de las tensiones de línea
    
    def matrix_vline(self):
        n=self.c;
        m=self.f;
        self.Vline=np.zeros((2*m*n,1))
        
#Método llamado por la clase para crear la matriz E en este caso, se debe pasar
#por parámetro un vector con las tensiones de linea aplicadas en el siguiente orden:
#           *Vappw1
#           *Vappw2
#           *Vappb1
#           *Vappb2                   
    
    def matrix_E(self, vapp):
        n=self.c;
        m=self.f;
        self.E=np.zeros((2*n*m,1))

#En primer lugar, se almacenan los valores correspondientes aplicadas en las
#Word Lines
 
        for i in range(m):
                    self.E[i*n,0]=vapp[i]/self.Rwl1;
                    self.E[(i*n)+(n-1),0]=vapp[i+m]/self.Rwl2;
                
                
#Luego se almacenan en la misma matriz los valores correspondientes a
#las tensiones aplicadas en las Bit Lines

        for i in range(n):
                    self.E[i*m+n*m,0]=-vapp[i+2*m]/self.Rbl1;
                    self.E[(i*m)+(m-1)+n*m,0]=-vapp[i+n+2*m]/self.Rbl2;
                    
                 
#Mapeo de Vline en dos matrices n*m
    
    def vline_mapping(self):
        n=self.c;
        m=self.f;
        
        self.VlineWL=np.zeros((m,n));
        self.VlineBL=np.zeros((m,n));
        
        for i in range (m):
            for j in range (n):
                self.VlineWL[i,j]=self.Vline[(i*n)+j];
                self.VlineBL[i,j]=self.Vline[(i*n)+j+n*m];
                


#Método utilizado por la clase para calcular la matriz inversa de ABCD
        
    def matrix_i_ABCD(self):
        self.ABCD_inv=np.linalg.inv(self.ABCD);          

#Método utilizado para la resolución directa
    #@lru_cache(maxsize=16)    
    def matrix_solve(self):
        self.Vline=sp.spsolve(self.ABCD,self.E);     

#Método utilizado para obtener las corrientes de salida
        
    def get_output(self):
        out=[];
        out=np.zeros(self.c);
        if(self.outputside=="down"):
            for i in range (self.c):
                out[i]=(self.VlineBL[self.f-1,i]-self.Vapp[2*self.f+self.c+i])/self.Rbl2
        else:
            for i in range (self.c):
                out[i]=(self.VlineBL[0,i]-self.Vapp[2*self.f+i])/self.Rbl1
        return out;
 
#Método utilizado para resolver nuevamente el Crossbar
  
    def cross_solve(self,vapp):
        self.Vapp=vapp;
        self.matrix_E(vapp);
        # t0=time.time();
        self.Vline=sp.spsolve(self.ABCD,self.E);

        self.vline_mapping(); 
        # t1=time.time()-t0;
        # print("\nSolved in: ",t1,"seconds\n")

#Método que permite calcular la matriz de factores de calibración
    def make_cal_factor(self):
        self.call=np.zeros((self.f,self.c));
        
        for i in range(self.f):
            for j in range(self.c):
                self.call[i,j]=self.Vapp[i]/self.VlineWL[i,j];

#Método de calibración de elementos Rij
                 
    def calibrate(self,Vcal,crit_val):
        
        it_num=1;
        
        max_cal_val=crit_val;
        min_cal_val=1;
        prom_cal_val=0;
        prom_cont=0;

        self.cross_solve(Vcal);
        self.cal=np.zeros((self.f,self.c));
        
        new_rij=0;
        
        for i in range(self.f):
            for j in range(self.c):
                self.cal[i,j]=self.Vapp[i]/self.VlineWL[i,j];
                if (self.cal[i,j]>0.001):
                    if(self.cal[i,j]>max_cal_val):
                        max_cal_val=self.cal[i,j];
                    if(self.cal[i,j]<min_cal_val):
                        min_cal_val=self.cal[i,j];
                    prom_cal_val=prom_cal_val+self.cal[i,j];
                    prom_cont=prom_cont+1;
                    new_rij=self.Rij[i,j]/self.cal[i,j];
                    if (new_rij>self.LRS and new_rij<self.HRS):
                        self.Rij[i,j]=new_rij;
       
        self.matrix_a();
        self.matrix_b();
        self.matrix_c();
        self.matrix_d();
        self.matrix_abcd();
        
        aux=[]
        aux=self.cal;
        criterion=1;
        
        if(prom_cont!=0):
            prom_cal_val=prom_cal_val/prom_cont;
        print("Loop de Cal. N°: ",it_num,"\nFactor de Cal. Máximo: ",max_cal_val,"\nFactor de Cal. Mínimo: ",min_cal_val,"\nFactor Cal. Promedio: ",prom_cal_val,"\n");
        
        self.cross_solve(Vcal);
        
        dif=[];
        past_dif=[];
        increasing_dif_counter=[];
        dif=np.zeros((self.f,self.c));
        past_dif=np.zeros((self.f,self.c));
        increasing_dif_counter=np.zeros((self.f,self.c));        
        new_rij=0;
       
        while (criterion==1):
            criterion=0;
            max_cal_val=crit_val;
            min_cal_val=1;
            prom_cal_val=0;
            prom_cont=0;
            dif_max=-999999;
            dif_prom=0;
            dif_min=9999;

            for i in range(self.f):
                for j in range(self.c):
                    past_dif[i,j]=dif[i,j];
                    dif[i,j]=np.absolute(aux[i,j]-(self.Vapp[i]/self.VlineWL[i,j]));
                    self.cal[i,j]=self.Vapp[i]/self.VlineWL[i,j];

                    if(dif[i,j]>dif_max):
                        dif_max=dif[i,j];
                    if(dif[i,j]<dif_min):
                        dif_min=dif[i,j];
                    dif_prom=dif_prom+dif[i,j];
                    if(self.cal[i,j]>max_cal_val):
                        max_cal_val=self.cal[i,j];
                    if(self.cal[i,j]<min_cal_val):
                        min_cal_val=self.cal[i,j];
                    prom_cal_val=prom_cal_val+self.cal[i,j];
                    prom_cont=prom_cont+1;

                    
                    if(dif[i,j]>crit_val):
                        if(dif[i,j]>past_dif[i,j]):                    
                            increasing_dif_counter[i,j]=increasing_dif_counter[i,j]+1;
                            if(increasing_dif_counter[i,j]==20):
                                return -1;
                        else:
                            increasing_dif_counter[i,j]=0;
                        
                        self.cal[i,j]=self.Vapp[i]/self.VlineWL[i,j];
                        if (self.cal[i,j]>0.001):
                            new_rij=self.Rij[i,j]/self.cal[i,j]
                            if (new_rij>self.LRS and new_rij<self.HRS):
                                self.Rij[i,j]=new_rij;
                        criterion=1;
 

            it_num=it_num+1;
            prom_cal_val=prom_cal_val/prom_cont;
            dif_prom=dif_prom/prom_cont;
            print("Loop de Cal. N°: ",it_num,"\nFactor de Cal. Máximo: ",max_cal_val,"\nFactor de Cal. Mínimo: ",min_cal_val,"\nFactor Cal. Promedio: ",prom_cal_val,"\n");
            print("Distancia máxima al valor previo: ",dif_max,"\nDistancia mínima al valor previo: ",dif_min,"\nPromedio de las diferencias: ",dif_prom,"\n")

            if(criterion==0):
                return it_num;
            

            
 
            
            self.matrix_a();
            self.matrix_b();
            self.matrix_c();
            self.matrix_d();
            self.matrix_abcd();
            self.cross_solve(Vcal);
            aux=self.cal;

        return it_num;
            
    
                     
              
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

# Rij=[];
# Vapp=[];

# Rij=np.zeros((f,c));

# for i in range (f):
#     for j in range (c):
#         if (i==j):
#             Rij[i,j]=100e3;
#         else:
#             Rij[i,j]=10e6;


# Vapp=np.zeros((2*(f+c)));
# # Vapp[3]=0.5;
# Vapp[8]=0.5;
# #Vapp[48]=0.5;

    
# mycrossbar = Crossbar(Rij,10,"SISO",Vapp);

# mycrossbar.make_call_factor();
# mycrossbar.callibrate();

# iout=[]

# iout=mycrossbar.get_output();


