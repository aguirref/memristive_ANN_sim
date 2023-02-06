*Memristor subcircuit developed by Chang et al.

* Connections:
* TE:  Top electrode
* BE:  Bottom electrode
* XSV: External connection to plot state variable
*      that is not used otherwise
.global gnd!
.SUBCKT memdiode p n H0=0

* Parameters:
* alpha:      Prefactor for Schottky barrier
* beta:       Exponent for Schottky barrier
* gamma:      Prefactor for tunneling
* delta:      Exponent for tunneling
* xmax:       Maximum value of state variable
* xmin:       Minimum value of state variable
* drift_bit:  Binary value to switch the ionic drift in (1)
*             or out (0) of the equation
* lambda:     State variable multiplier
* eta1, eta2: State variable exponential rates
* tau:        Diffusion coefficient

.param alpha=0.1e-5 
.param beta=0.5 
.param gamma=8e-5
.param delta=1.7 
.param xmax=1 
.param xmin=0
.param drift_bit = 0.4 
.param lambda=1 
.param eta1=0.006 
.param eta2=12.5 
.param tau=3e2
.param cp=1

* Auxiliary functions to limit the range of x
.param sign2(var) = '(sgn(var)+1)/2'
.param trun(var1,var2) = 'sign2(var1)*sign2(xmax-var2)+sign2(-var1)*sign2(var2-xmin)'

* Memristor IV Relationship
Gm p n cur='(1-cp*V(H))*alpha*(1-exp(-beta*V(p,n)))+(cp*V(H))*gamma*sinh(delta*V(p,n))'

* Rate equation for state variable
Gx gnd! H cur='trun(V(p,n),cp*V(H))*lambda*(eta1*sinh(eta2*V(p,n))-drift_bit*cp*V(H)/tau)'
Cpvar H gnd! C='cp' IC='H0'

.ENDS memdiode
