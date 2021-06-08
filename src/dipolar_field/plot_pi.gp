# set term epslatex color colortext size 6.4,8.0
# set output "dipolo_0.tex"

set term qt size 800, 800


set multiplot layout 3,2 title "\\textbf{Distribuzioni di probabilità per} $\\mathbf{\\theta = 0}$"
unset key 
set grid

set xlabel "$z^*$"
set ylabel "$\\pi(z^*)$"
set title '$E^* = 0.1$'
plot 'pi.csv' using 1:2 w lines lc 7 lt 7

set xlabel "$z^*$"
set ylabel "$\\pi(z^*)$"
set title '$E^* = 0.3$'
plot 'pi.csv' using 1:3 w lines lc 2 lt 7

set xlabel "$z^*$"
set ylabel "$\\pi(z^*)$"
set title '$E^* = 0.5$'
plot 'pi.csv' using 1:4 w lines lc 4 lt 7

set xlabel "$z^*$"
set ylabel "$\\pi(z^*)$"
set title '$E^* = 0.1$'
plot 'pi.csv' using 1:5 w lines lc 9 lt 7

set xlabel "$z^*$"
set ylabel "$\\pi(z^*)$"
set title '$E^* = 2.0$'
plot 'pi.csv' using 1:6 w lines lc 11 lt 7

set xlabel "$z^*$"
set ylabel "$\\pi(z^*)$"
set title '$E^* = 5.0$'
plot 'pi.csv' using 1:7 w lines lc 24 lt 7

pause -1