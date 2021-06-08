set datafile separator '\t'

set title "Quantum 3D scattering"
set xlabel 'E [meV]'
set ylabel '{/Symbol S}_{tot} [A^2]' offset 1,0 rotate by 0
unset key
set grid 
plot 'sigma_tot.csv' using 1:2 lt 7 with lines 

pause -1