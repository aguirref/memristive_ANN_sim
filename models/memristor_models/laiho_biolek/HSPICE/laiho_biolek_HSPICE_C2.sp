* SPICE model for equations proposed by Dr. Mika Laiho et al.
* with Biolek window for memristor boundaries
* Connections:
* TE: Top electrode
* BE: Bottom electrode
* XSV: External connection to plot state variable
* that is not used otherwise
***************
.global gnd!

.subckt memdiode p n H0=1e-6

.param a1=9e-5 
.param b1=1.5 
.param a2=0.7e-4 
.param b2=1.8 
.param c1=4.5e-2 
.param d1=9 
.param c2=7.5e-1 
.param d2=1.8 
.param p=1
.param CH0=1

* Hyperbolic sine IV relationship
.param IVRel(V1,V2) = 'V1 >= 0 ? a1*V2*sinh(b1*V1) : a2*V2*sinh(b2*V1)'
* unitary step function
.param STP(V) = '(sgn(V)+1)/2'
* Equation for state variable
.param SV(V1) = 'V1 >= 0 ? c1*sinh(d1*V1) : c2*sinh(d2*V1)'
* Biolek window function
.param f(V1,I1) = '1-pow((V1-STP(-I1)),(2*p))'

* Current source representing memristor

Gm p n cur='IVRel(V(p,n),V(H,gnd!))'

* Circuit to determine value of state variable
Gx   gnd! H    cur='SV(V(p,n))*f(V(H,gnd!),I(Gm))'
Cx   H    gnd! C='CH0' IC='H0'

.ends memdiode
