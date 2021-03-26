set datafile separator "\t"

set xrange [-5:5]
set yrange [-5:5]
set zrange [-5:5]

splot "Output/trajectory.csv" using 2:3:4