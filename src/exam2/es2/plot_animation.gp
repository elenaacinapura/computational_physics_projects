reset
set term qt 0
set xrange [-1.5:1.5]
set yrange [-0.5:2.5]
set grid

set title 'Densità di probabilità'

col = 100
do for [ii=4:col + 2] {
    plot 'animation.csv' using 1:ii w l ls 1 notitle ,\
        'animation.csv' using 1:2 w l ls 7 notitle
    pause 0.001
}