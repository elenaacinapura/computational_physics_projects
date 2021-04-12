set term pngcairo

set output "T.png"

set datafile separator "\t"

unset key

plot "T.csv" using 1:2 lc 7 w lines