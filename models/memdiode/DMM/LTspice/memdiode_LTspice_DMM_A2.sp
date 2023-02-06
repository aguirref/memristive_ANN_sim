***************
.subckt memdiode p n H0=0
.param beta=5.000000000000e-01
.param EI=0
.param T0s=1.202042000000e+02
.param V0s=9.789000000000e-02
.param T0r=NaN
.param V0r=NaN
.param imax=9.576000000000e-05
.param alphamax=9.999000000000e-01
.param rsmax=3.000000000000e+01
.param imin=1.800000000000e-06
.param alphamin=1.000000000000e+00
.param rsmin=3.000000000000e+01
*Auxiliary functions
.param I0(x)='imax*x+imin*(1-x)'
.param A(x)='alphamax*x+alphamin*(1-x)'
.param Rss(x)='rsmax*x+rsmin*(1-x)'

*H-V
EV A 0 vol=1
RH H A R='T0s*exp(-V(p,n)/V0s)'
RD H 0 R='T0r*exp(V(p,n)/V0r)'
CH H 0 1 IC='H0'

*I-V
 RS p D R='Rss(V(H))'
GD D n cur='I0(V(H))*(exp(beta*A(V(H))*V(D,n))-exp(-(1-beta)*A(V(H))*V(D,n)))+EI'

.ends memdiode
