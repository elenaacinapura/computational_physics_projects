set term epslatex size 5.0, 9.0 color colortext solid font ",7"
set output "levels.tex"

# set term qt size 650, 1000

set title "Stati legati in approssimazione Bohr-Sommerfeld"

set multiplot layout 3,1 scale 0.8,1 title "Stati legati in approssimazione Bohr-Sommerfeld" font ",17"


set key right center
set xrange [-3:3]
set xlabel '$x^*$'
set ylabel "$E^*$" offset 1,0 rotate by 0
set ytic 0.2


xi = 0.025
array E [4] = [0.787, 0.431, 0.172, 0.02]
set title sprintf("$\\xi = $ %.03f", xi) font ",14"
set for [i=1:4] arrow i from -acosh((1.0/E[i])**(1.0/4)),-E[i] to acosh((1/E[i])**(1.0/4)),-E[i] nohead lw 2 lc 7
plot -cosh(x)**(-4) w lines lc 0 title "$V(x)$"

xi = 0.005
array E [9] = [0.902, 0.720, 0.556, 0.411, 0.285, 0.180, 0.097, 0.037, 0.004]
set title sprintf("$\\xi = $ %.03f", xi)
set for [i=1:9] arrow i from -acosh((1.0/E[i])**(1.0/4)),-E[i] to acosh((1/E[i])**(1.0/4)),-E[i] nohead lw 2 lc 7
plot -cosh(x)**(-4) w lines lc 0 title "$V(x)$"

xi = 0.00025
array E [13] = [0.930, 0.797, 0.674, 0.559, 0.455, 0.359, 0.274, 0.199, 0.135, 0.082, 0.0409, 0.0133, 0.0005]
set title sprintf("$\\xi = $ %.05f", xi)
set for [i=1:13] arrow i from -acosh((1.0/E[i])**(1.0/4)),-E[i] to acosh((1/E[i])**(1.0/4)),-E[i] nohead lw 2 lc 7
plot -cosh(x)**(-4) w lines lc 0 title "$V(x)$"


unset multiplot