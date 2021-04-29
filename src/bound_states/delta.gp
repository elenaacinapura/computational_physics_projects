set datafile separator '\t'
set term qt 0
set xlabel "E^*"
set ylabel "{/Symbol D}" offset 1,0 rotate by 0


set arrow from -1,0 to 0,0 nohead
set title 'Cosh potential'
plot 'delta_cosh.csv' using 1:2 w lines notitle

set term qt 1
set title 'Lennard Jones potential'
set yrange [-20:20]
plot 'delta_lj.csv' using 1:2 w lines notitle
pause -1

