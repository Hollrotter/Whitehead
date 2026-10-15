set terminal png size 400, 300

# set format y "%.0e"
# set logscale y 10
set xlabel '# Knoten'

set grid
set output '../png/SquareWingTime.png'

set multiplot layout 2,1 rowsfirst
set ylabel 'cL'
set yrange [0.44:0.6]
set ytics 0.44,0.04
set lmargin 10
plot '../Data/SquareWing' u 1:5 notitle with lines lw 2 lc rgb 'red'
set ylabel 't/s'
set yrange [0:160]
set ytics 0,40
set lmargin 10
plot '../Data/SquareWing' u 1:7 notitle with lines lw 2 lc rgb 'red'
unset multiplot