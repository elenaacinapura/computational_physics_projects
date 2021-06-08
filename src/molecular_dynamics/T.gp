set term qt

set datafile separator "\t"

unset key

set title "Andamento della temperatura"
set xlabel 't'
set ylabel 'T' offset 1,0 rotate by 0
plot "T.csv" using 1:2 lc 7 w lines

pause -1