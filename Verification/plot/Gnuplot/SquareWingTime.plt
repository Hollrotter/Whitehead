set terminal png size 400, 300

# set format y "%.0e"
# set logscale y 10
set xlabel '# Knoten'
set ylabel 't/s'
set xrange [0:16]
set yrange [0:5]
set ytics 0,1
set lmargin 10

set grid
set output '../png/SquareWingTime.png'

plot '../Data/SquareWing' u 1:7 notitle with lines lw 2 lc rgb 'red'