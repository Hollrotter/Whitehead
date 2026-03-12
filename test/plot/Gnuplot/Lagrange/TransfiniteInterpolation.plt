set terminal png size 1000, 1000

set output '../../png/Lagrange/TransfiniteInterpolation.png'
set palette defined (0 "black", 0.337 "#000080", 0.675 "#0080ff", 1 "white")
set palette maxcolor 100
set xlabel 'x'
set ylabel 'y'

# unset border
# unset xtics
# unset ytics

# set view map
set pm3d map
set pm3d interpolate 10,10 corners2color mean

# splot '../../Data/Lagrange/TransfiniteInterpolation' notitle with lines lw 2 lc rgb 'black'
splot '../../Data/Lagrange/TransfiniteInterpolation' u 1:2:6 notitle