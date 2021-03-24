set datafile separator "\t"

unset key

plot "Output/energies.csv" using 1:2 lc 7 w lines