set terminal png size 350, 300

# set format y "%.0e"
# set logscale y 10
set xlabel '# Knoten'
set ylabel 'cL'
set grid
set output '../png/SquareWing.png'

plot '../Data/SquareWing' u 1:5 notitle with lines lw 2 lc rgb 'red'