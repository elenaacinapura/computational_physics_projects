reset

set term qt 0

####### SPECTRUM ##########
set xrange[0:0.5]
plot "spectrum.csv" using 1:2 w l lt 7
pause -1

####### PROBABILITY ##########
set xrange [0:2048]
#set yrange [0:0.5]
plot "probability.csv" using 1:2 w l lt 7
pause -1
