.global gnd! vss! vdd!

.subckt ideal_comparator p n out 

.param gain=1e8

EV innver_voltage gnd! vol='(V(p)-V(n)*gain'

gi1 out vdd! vcr pwl(1) inner_voltage gnd! 1e-3,1000Meg 0,1
gi2 out vss! vcr pwl(1) gnd! inner_voltage 1e-3,1000Meg 0,1

.ends ideal_comparator
