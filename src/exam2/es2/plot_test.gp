reset

set term qt 0

####### SPECTRUM ##########
set xrange[0:1.5]
plot "spectrum.csv" using 1:2 w l lt 7
pause -1

####### PROBABILITY ##########
#set xrange [0:3000]
#set yrange [0:1]
#plot "probability.csv" using 1:2 w l lt 7
#pause -1

######## ANIMATION ############
#set xrange [-1.5:1.5]
#set yrange [-0.5:2.5]
#set grid

#col = 100
#do for [ii=4:col + 2] {
#    plot 'animation.csv' using 1:ii w l ls 1 ,\
#        'animation.csv' using 1:2 w l ls 7
#    pause 0.001
#}