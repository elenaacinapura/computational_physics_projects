set datafile separator "\t"

plot "Output/energies.csv" using 1:4 title 'Total'      w points lc 7 ,\
#     "Output/energies.csv" using 1:2 title 'Kinetic'   w lines lc 4 ,\
#     "Output/energies.csv" using 1:3 title 'Potential' w lines lc 2 