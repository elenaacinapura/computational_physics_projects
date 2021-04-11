set term pngcairo

set output "energies.png"

set datafile separator "\t"

set title "Energie"

set xlabel "t*"
set ylabel "E*"

set arrow from 80, graph(0,0) to 80,graph(1,1) nohead

set label "Fine dell'equilibrazione" at 78,-800 rotate by 90

plot "energies.csv" using 1:4 title 'Totale'      w lines lc 6 lw 2,\
     "energies.csv" using 1:2 title 'Cinetica'    w lines lc 39 ,\
     "energies.csv" using 1:3 title 'Potenziale'  w lines lc 3,\

