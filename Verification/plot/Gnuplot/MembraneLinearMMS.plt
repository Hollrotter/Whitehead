set terminal png size 350, 300

set format y "%.0e"
set logscale y 10
set xlabel '# Knoten'
set ylabel 'Fehler'
set grid
set output '../png/MembraneLinearMMS.png'

plot '../Data/MembraneLinearMMS' u 1:3 notitle with lines lw 2 lc rgb 'red'