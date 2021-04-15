set datafile separator '\t'

set xlabel 'E [meV]'
set ylabel '{/Symbol S}_{tot}' offset 1,0 rotate by 0
unset key
set grid 
plot 'sigma_tot.csv' using 1:2 with lines lt 7 