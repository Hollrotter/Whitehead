set terminal png size 350, 300

set format y "%.0e"
set logscale y 10
set xlabel 'Number of nodes'
set ylabel 'RMS error'
set grid
set output '../png/MembraneNonlinearMMS.png'

plot '../Data/MembraneNonlinearMMS' u 1:2 title "x-direction" with lines lw 2 lc rgb 'green',\
     '../Data/MembraneNonlinearMMS' u 1:3 title "y-direction" with lines lw 2 lc rgb 'blue',\
     '../Data/MembraneNonlinearMMS' u 1:4 title "z-direction" with lines lw 2 lc rgb 'red'