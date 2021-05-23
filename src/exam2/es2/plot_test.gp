reset

set term qt 0

set xrange [-2:2]
set yrange [-0.5:2]


col = 100
do for [ii=4:col + 2] {
    plot 'animation.csv' using 1:ii w l ls 1 ,\
        'animation.csv' using 1:2 w l ls 7
    pause 0.01
}

#plot "test.csv" using 1:2 w l lt 3 title 'potential' ,\
#    "test.csv" using 1:3 w l lt 5 title 'initial' ,\
#    "test.csv" using 1:4 w l lt 7 title 'final'
#pause -1