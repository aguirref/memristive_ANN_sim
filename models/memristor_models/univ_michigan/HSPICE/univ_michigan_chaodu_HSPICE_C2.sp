****************** LTspice code for metal oxide memristors*****************
*Parameters:
*alpha is prefactor for Schottky barrier
*beta is exponent prefactor for Schottky barrier
*gamma is prefactor for tunneling
*delta is exponent prefactor for tunneling
***************************************************************************
.global gnd!
.SUBCKT memdiode p n H0=0

.param alpha=2e-7
.param beta=2.1 
.param gamma=1e-4 
.param delta=1.35 
.param wmax=1 
.param wmin=0

*State variable:
.param lambda=1e-2 
.param rhoc=12.5 
.param rhom=0.1
.param taul=1e3 
.param taus=1e3 
.param epsilon=0.01 
.param sigma=25
.param cc=1
.param cm=1

Cpvar1 c gnd! C='cc' IC='H0'
Cpvar2 m gnd! C='cm' IC=0.001

*rate equation considering the diffusion effect
Gc gnd! c cur='trunc1(V(p,n),cc*V(c))*(lambda*exp(epsilon*cc*V(c))*sinh(rhoc*V(p,n)))-(cc*V(c)-0.001)*(1/taul+sigma*cm*V(m)/taus)'
Gm gnd! m cur='trunc2(V(p,n),cm*V(m))*(lambda*sinh(rhom*abs(V(p,n))))-(cm*V(m)-0.001)*(cm*V(m)/taus)'

**********************************************************
*auxiliary functions to limit the range of w
.param sign2(var)='(sgn(var)+1)/2'
.param trunc1(var1,var2)='sign2(var1)*sign2(wmax-var2)*(1-exp(-(wmax-var2)/0.0001))+sign2(-var1)*sign2(var2-wmin)*(1-exp(-(var2-wmin)/0.0001))'
.param trunc2(var1,var2)='sign2(var1)*sign2(wmax-var2)*(1-exp(-(wmax-var2)/0.0001))+sign2(-var1)*sign2(var2-wmin)*(1-exp(-(var2-wmin)/0.0001))'

***************************************************************************
*Output:
Gw p n cur='(1-cc*V(c))*alpha*(1-exp(-beta*V(p,n)))+(cc*V(c))*gamma*sinh(delta*V(p,n))'

.ENDS memdiode
