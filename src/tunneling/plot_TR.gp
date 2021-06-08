# set term epslatex color colortext size 7.0,3.5
# set output "tunneling.tex"

set terminal qt size 1000, 500
set datafile separator "\t"

# Coefficiente di trasmissione

set multiplot layout 1,2 title "\\big{\\textbf{Effetto tunnel quantistico}}"
set grid

set xlabel "$E^*$"
set ylabel "$|T|^2$"
set key right bottom

set title "Coefficiente di trasmissione"
plot "T.csv" using 1:2 w lines lc 7 lw 2 title "$\\xi=0.01$" ,\
     "T.csv" using 1:4 w lines lc 11 lw 2 title "$\\xi=0.025$" ,\
     "T.csv" using 1:6 w lines lc 20 lw 2 title "$\\xi=0.005$"



# Coefficiente di riflessione 

set xlabel "$E^*$"
set ylabel "$|R|^2$"
set key right top

set title "Coefficiente di riflessione"
plot "T.csv" using 1:3 w lines lc 7 lw 2 title "$\\xi=0.01$" ,\
     "T.csv" using 1:5 w lines lc 11 lw 2 title "$\\xi=0.025$" ,\
     "T.csv" using 1:7 w lines lc 20 lw 2 title "$\\xi=0.005$"  

pause -1