# set terminal qt
set term epslatex color colortext 
set output "variational.tex"

set datafile separator "\t"

set xlabel "$\\alpha*"
set ylabel "E*" offset 2,1 rotate by 0

set object circle at 0.522,-0.797 size 0.01 fillcolor rgb "red" fillstyle solid
set object circle at 0.330,-0.904 size 0.01 fillcolor rgb "cyan" fillstyle solid
set object circle at 0.274,-0.931 size 0.01 fillcolor rgb "orange" fillstyle solid

set title "$E^*_{\\alpha^*}$ in funzione di $\\alpha^*$"
plot "E.csv" using 1:2 w lines lc 7  lw 2 title "$\\xi=0.025$" ,\
     "E.csv" using 1:3 w lines lc 11 lw 2 title "$\\xi=0.005$" ,\
     "E.csv" using 1:4 w lines lc 20 lw 2 title "$\\xi=0.0025$"