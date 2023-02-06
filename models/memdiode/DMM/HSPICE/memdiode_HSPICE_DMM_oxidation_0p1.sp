***************
.SUBCKT memdiode p n H0=0 x0=0.1

.param imax_slope=-0.07005643624874254
.param imax_intercept=-2.31749382490389
*.param imax=1
.param imax='pwr(10,linearmodel(x0,imax_slope,imax_intercept))'

.param imin_slope=-2.233493992404483
.param imin_intercept=-4.186285859587773
*.param imin=1
.param imin='pwr(10,linearmodel(x0,imin_slope,imin_intercept))'

.param alphamax_slope=1.1384749367818683
.param alphamax_intercept=2.430231899198936
*.param alphamax=1
.param alphamax='linearmodel(x0,alphamax_slope,alphamax_intercept)'

.param alphamin_slope=2.2440796366498863
.param alphamin_intercept=2.723028742076817
*.param alphamin=1
.param alphamin='linearmodel(x0,alphamin_slope,alphamin_intercept)'

.param vset_slope=1.0607990611668159
.param vset_intercept=0.791481114033967
*.param vset=1
.param vset='linearmodel(x0,vset_slope,vset_intercept)'

.param vres_slope=0.8000000000000002
.param vres_intercept=-1.6499999999999997
*.param vres=1
.param vres='linearmodel(x0, vres_slope, vres_intercept)'

.param etaset_a=1.5000000000503808
.param etaset_b=1.6499999999145025
.param etaset_c=2.449999999914503
*.param etaset=1
.param etaset='pwr(10,logmodel(x0,etaset_a,etaset_b,etaset_c))'

.param etares_slope=0.48077237638582676
.param etares_intercept=0.2798181704934635
*.param etares=1
.param etares='pwr(10,linearmodel(x0, etares_slope, etares_intercept))'

.param EI=1e-10
.param isb=1
.param beta=0.5
.param ch0=1
.param rsmax=1
.param rsmin=1
.param vt=0.1
.param gam0=0
.param gam=0

*Auxiliary functions
.param I0(x)='imin+(imax-imin)*x'
.param A(x)='alphamin+(alphamax-alphamin)*x'
.param RS(x)='rsmin+(rsmax-rsmin)*x'
.param linearmodel(x, slope, intercept)='slope*x+intercept'
.param logmodel(x, a, b, c)='a*log10(x*b)+c'
*.param TS(x)='exp(-etaset*(x-VSNAPB(I(RS))))'
.param TS(x)='exp(-etaset*(x-vset))'
.param TR(x)='exp(etares*ISF(V(H))*(x-vres))'
.param VSNAPB(x)='x>isb ? vt : vset'
.param ISF(x)='gam==0 ? 1 : (pwr(x,gam)-gam0)'

*H-V
GI gnd! H cur='V(p,n)>=0 ? (1-V(H))/TS(V(D,n)) : -V(H)/TR(V(D,n))'
CH H gnd! C=1 ic='H0'
*RH  H gnd! R=1E9

*EV A gnd! vol='V(p,n)>=0 ? 1 : 0'
*RH A H R='V(p,n)>=0 ? TS(V(D,n)) : TR(V(D,n))'
*CH H gnd! C=1 ic='H0'

*I-V
RS p D r='RS(V(H))'
GD D n cur='I0(V(H))*(exp(beta*A(V(H))*V(D,n))-exp(-(1-beta)*A(V(H))*V(D,n)))+EI'

.ENDS
