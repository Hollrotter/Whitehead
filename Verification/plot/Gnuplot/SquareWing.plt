set terminal png size 350, 300

# set format y "%.0e"
# set logscale y 10
set xlabel '# Knoten'

set grid
set output '../png/SquareWing.png'

set multiplot layout 2,1 rowsfirst
set ylabel 'cL'
set yrange [0.08:0.1]
set ytics 0.08,0.005
plot '../Data/SquareWing' u 1:5 notitle with lines lw 2 lc rgb 'red'
set ylabel 't/s'
set yrange [0:100]
set ytics 0,20
plot '../Data/SquareWing' u 1:7 notitle with lines lw 2 lc rgb 'red'
unset multiplot