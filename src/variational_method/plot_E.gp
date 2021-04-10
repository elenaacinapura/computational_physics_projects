set terminal pngcairo 
set datafile separator "\t"

# Energie in funzione di alpha

set output "E.png"

set xlabel "alpha*"
set ylabel "E*"

set title "Energie in funzione di alpha"
plot "E.csv" using 1:2 w lines lc 7  title "xi=0.2" ,\
     "E.csv" using 1:3 w lines lc 11 title "xi=0.05" ,\
     "E.csv" using 1:4 w lines lc 20  title "xi=0.01"