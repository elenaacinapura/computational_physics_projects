set term pngcairo

set output "g.png"

set datafile separator "\t"

unset key

set title "Funzione di distribuzione a coppie"

set xlabel "r*"
set ylabel "g(r)"

plot "g.csv" using 1:2 lc 14 lt 7 w linespoint