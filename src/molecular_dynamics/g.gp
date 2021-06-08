set term epslatex color colortext size 6.0,7.0
set output "g.tex"
# set term qt size 900,1000

set datafile separator "\t"

unset key

set multiplot layout 3,2 title "\\textbf{Funzione di distribuzione a coppie}"

set xlabel "$r^*$"
set ylabel "$g(r)$" offset 1,0 rotate by 0

set title "$\\rho^* = 0.5 \\quad T^* = 0.5$"
plot "g_0.csv" using 1:2 lc 14 lt 7 w linespoint

set title "$\\rho^* = 0.5 \\quad T^* = 0.8$"
plot "g_1.csv" using 1:2 lc 14 lt 7 w linespoint

set title "$\\rho^* = 0.6 \\quad T^* = 0.5$"
plot "g_2.csv" using 1:2 lc 14 lt 7 w linespoint

set title "$\\rho^* = 0.6 \\quad T^* = 0.8$"
plot "g_3.csv" using 1:2 lc 14 lt 7 w linespoint

set title "$\\rho^* = 0.7 \\quad T^* = 0.5$"
plot "g_4.csv" using 1:2 lc 14 lt 7 w linespoint

set title "$\\rho^* = 0.7 \\quad T^* = 0.8$"
plot "g_5.csv" using 1:2 lc 14 lt 7 w linespoint