set terminal png size 650, 450

# set format y "%.0e"
# set logscale y 10
set xlabel 'Anzahl Knoten (Strukturmodell)'
set ylabel 't / s'
set xrange [6:16]
set yrange [0:5]
set ytics 0,1
set lmargin 10

set grid
set output '../png/SquareWing_t.png'

plot '../Data/SquareWing' u 1:7 notitle with lines lw 2 lc rgb 'blue'