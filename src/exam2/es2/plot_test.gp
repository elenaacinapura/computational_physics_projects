reset

set term qt 0

set xrange [-2:2]
set yrange [-0.5:2]

plot "ground.csv" using 1:2 w l lt 3 ,\
    "ground.csv" using 1:3 w l lt 5
pause -1