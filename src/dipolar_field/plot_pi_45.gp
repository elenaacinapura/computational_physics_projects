set term pngcairo size 1800,1000 font "Helvetica,20"

set output "pi_45.png" 

set multiplot layout 2,3 scale 1,1

set title "Theta = pi/4"

set xlabel "z"
set ylabel "pi(z)"
plot 'pi_45.csv' using 1:2 w lines lc 7 lt 7 title 'E = 0.1'
set xlabel "z"
set ylabel "pi(z)"
plot 'pi_45.csv' using 1:3 w lines lc 2 lt 7 title 'E = 0.3'
set xlabel "z"
set ylabel "pi(z)"
plot 'pi_45.csv' using 1:4 w lines lc 4 lt 7 title 'E = 0.5'
set xlabel "z"
set ylabel "pi(z)"
plot 'pi_45.csv' using 1:5 w lines lc 9 lt 7 title 'E = 1.0'
set xlabel "z"
set ylabel "pi(z)"
plot 'pi_45.csv' using 1:6 w lines lc 11 lt 7 title 'E = 2.0'
set xlabel "z"
set ylabel "pi(z)"
plot 'pi_45.csv' using 1:7 w lines lc 17 lt 7 title 'E = 5.0'
