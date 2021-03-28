set datafile separator "\t"

set arrow from 40, graph(0,0) to 40,graph(1,1) nohead
set label "End of equilibration" at 38,-800 rotate by 90

plot "Output/energies.csv" using 1:4 title 'Total'      w lines lc 7 ,\
     "Output/energies.csv" using 1:2 title 'Kinetic'    w lines lc 4 ,\
     "Output/energies.csv" using 1:3 title 'Potential'  w lines lc 2 ,\

