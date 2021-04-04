set datafile separator "\t"

set xrange [-2.7:2.7]
set yrange [-2.7:2.7]
set zrange [-2.7:2.7]

splot "Output/traj.csv" using 2:3:4
pause 0.2
reread 