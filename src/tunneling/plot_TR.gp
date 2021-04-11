set terminal pngcairo 
set datafile separator "\t"

# Coefficiente di trasmissione

set output "T.png"

set xlabel "E*"
set ylabel "|T|^2"

set title "Coefficiente di trasmissione"
plot "T.csv" using 1:2 w lines lc 7  title "xi=0.2" ,\
     "T.csv" using 1:4 w lines lc 11 title "xi=0.05" ,\
     "T.csv" using 1:6 w lines lc 20  title "xi=0.01"



# Coefficiente di riflessione 

set output "R.png"

set xlabel "E*"
set ylabel "|R|^2"

set title "Coefficiente di riflessione"
plot "T.csv" using 1:3 w lines lc 7  title "xi=0.2" ,\
     "T.csv" using 1:5 w lines lc 11 title "xi=0.05" ,\
     "T.csv" using 1:7 w lines lc 20  title "xi=0.01"  