set term pngcairo

set output "trajectory.png"

unset key

plot "trajectory.csv" using 1:2 w lines, "<echo '0 0'" w points lt 7, "<echo '0 1'" w points lt 7