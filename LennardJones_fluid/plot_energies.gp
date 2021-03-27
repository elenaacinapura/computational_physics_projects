set datafile separator "\t"

plot "Output/energies.csv" using 1:2 w lines title 'Kinetic',\
     "Output/energies.csv" using 1:3 w lines title 'Potential',\
     "Output/energies.csv" using 1:4 w lines title 'Total'