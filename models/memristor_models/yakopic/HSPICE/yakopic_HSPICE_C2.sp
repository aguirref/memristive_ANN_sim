* SPICE model for memristive devices
* Created by Chris Yakopcic
* Last Update: 12/21/2011
*
* Connections:
* TE - top electrode
* BE - bottom electrode
* XSV - External connection to plot state variable
* that is not used otherwise
***************
.global gnd!

.subckt memdiode p n HO=1e-6
* Fitting parameters to model different devices
* a1, a2, b: Parameters for IV relationship
* Vp, Vn: Pos. and neg. voltage thresholds
* Ap, An: Multiplier for SV motion intensity
* xp, xn: Points where SV motion is reduced
* alphap, alphan: Rate at which SV motion decays
* xo: Initial value of SV
* eta: SV direction relative to voltage

.param a1=7e-5 
.param a2=7e-5 
.param b=1.9 
.param Vp=0.4 
.param Vn=0.3
.param Ap=1e-0 
.param An=1e-0 
.param xp=0.9 
.param xn=0.7
.param alphap=3
.param alphan=8
.param eta=1
.param CH0=1

* Multiplicative functions to ensure zero state
* variable motion at memristor boundaries
.param wp(V) = '(xp-V)/(1-xp)+1'
.param wn(V) = 'V/(1-xn)'
* Function G(V(t)) - Describes the device threshold
.param G(V) = 'V <= Vp ? (V >= -Vn ? 0 : -An*(exp(-V)-exp(Vn))): Ap*(exp(V)-exp(Vp))'
* Function F(V(t),x(t)) - Describes the SV motion
.param F(V1,V2) = 'eta*V1 >= 0 ? (V2 >= xp ? (exp(-alphap*(V2-xp))*wp(V2)) : 1) : (V2 <= (1-xn) ? exp(alphan*(V2+xn-1))*wn(V2) : 1)'
* IV Response - Hyperbolic sine due to MIM structure
.param IVRel(V1,V2) = 'V1 >= 0 ? a1*V2*sinh(b*V1) : a2*V2*sinh(b*V1)'
* Circuit to determine state variable
* dx/dt = F(V(t),x(t))*G(V(t))

Cx h gnd! C='CH0' IC='H0'
Gx gnd! h cur='eta*F(V(p,n),V(h,gnd!))*G(V(p,n))'

* Current source for memristor IV response

Gm p n cur='IVRel(V(p,n),V(h,gnd!))'

.ends memdiode
