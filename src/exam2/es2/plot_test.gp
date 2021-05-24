reset
set term qt 0

####### PROBABILITY ##########
set xrange [0:2048]
set title 'Probability p(t)'
set grid
set xlabel 'time'
set ylabel 'probability'
plot "probability.csv" using 1:2 w l lt 7 notitle
pause -1

####### SPECTRUM ##########
reset
set title 'Spectrum of p(t)'
set xrange[-0.5:0.5]
plot "spectrum.csv" using 1:2 w l lt 7 notitle
pause -1

