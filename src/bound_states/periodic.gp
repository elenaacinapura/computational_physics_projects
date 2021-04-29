set datafile separator '\t'

set term qt
set title 'Periodic square potential'
set xlabel "K"
set ylabel "E_n" offset 1,0 rotate by 0
plot 'periodic.csv' using 1:2 w lines title 'n = 0','periodic.csv' using 1:3 w lines title 'n = 1', 'periodic.csv' using 1:4 w lines title 'n = 2'
pause -1