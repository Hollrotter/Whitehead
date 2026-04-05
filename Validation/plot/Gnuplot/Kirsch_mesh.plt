set terminal svg size 1375, 500

set output '../svg/Kirsch_mesh.svg'

# set xlabel 'x/m'
# set ylabel 'y/m'

unset border
unset xtics
unset ytics

set view map

splot '../Data/Kirsch_s1_0' u 1:2:(0*$3) notitle with lines lw 2 lc rgb 'black',\
      '../Data/Kirsch_s1_1' u 1:2:(0*$3) notitle with lines lw 2 lc rgb 'black',\
      '../Data/Kirsch_s1_2' u 1:2:(0*$3) notitle with lines lw 2 lc rgb 'black',\
      '../Data/Kirsch_s1_3' u 1:2:(0*$3) notitle with lines lw 2 lc rgb 'black',\
      '../Data/Kirsch_s1_4' u 1:2:(0*$3) notitle with lines lw 2 lc rgb 'black'