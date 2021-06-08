set term epslatex color colortext
set output "energies.tex"
# set term qt

set datafile separator "\t"

set title "\\textbf{Andamento temporale delle energie}"
set key right bottom

set xlabel "$t^*$"
set ylabel "$E^*$" offset 1,0 rotate by 0

set label "Fine equilibrazione" font ",5" at 78,10 rotate by 90
set arrow from 80, graph(0,0) to 80,graph(1,1) nohead


plot "energies.csv" using 1:4 title 'Totale'      w lines lc 39 lw 2,\
     "energies.csv" using 1:2 title 'Cinetica'    w lines lc 6 ,\
     "energies.csv" using 1:3 title 'Potenziale'  w lines lc 3,\

pause -1
