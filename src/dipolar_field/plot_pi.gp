set term pngcairo

set output "pi.png"

unset key

plot 'pi.csv' using 1:2 w lines lc 7 lt 7
