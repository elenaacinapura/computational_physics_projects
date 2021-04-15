set datafile separator '\t'

set xlabel "E^*"
set ylabel "{/Symbol D}" offset 1,0 rotate by 0

set arrow from -1,0 to 0,0 nohead
plot 'delta.csv' using 1:2 w lines notitle