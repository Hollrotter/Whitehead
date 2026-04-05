set terminal svg size 1000, 1000

set output '../svg/CircularMembrane_mesh.svg'

# set xlabel 'x/m'
# set ylabel 'y/m'

unset border
unset xtics
unset ytics

set view map

splot '../Data/CircularMembrane_z_0' u 1:2:(0*$3) notitle with lines lw 2 lc rgb 'black',\
      '../Data/CircularMembrane_z_1' u 1:2:(0*$3) notitle with lines lw 2 lc rgb 'black',\
      '../Data/CircularMembrane_z_2' u 1:2:(0*$3) notitle with lines lw 2 lc rgb 'black'