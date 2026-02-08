set terminal png size 1000, 1000

set output '../../png/Structure/vy.png'
# set cbrange [0:0.0075]
# set palette defined (0 "black", 0.01 "#000080", 0.25 "#0080ff", 1 "white")
set palette defined (0  0.0 0.0 1.0, \
                     1  0.0 0.698 1.0, \
                     2  0.0 1.0 1.0, \
                     3  0.0 1.0 0.698, \
                     4  0.0 1.0 0.0, \
                     5  0.698 1.0 0.0, \
                     6  1.0 1.0 0.0, \
                     7  1.0 0.698 0.0, \
                     8  1.0 0.0 0.0 )
set palette maxcolor 100

set xlabel 'x/m'
set ylabel 'y/m'

set view map
set pm3d interpolate 10,10 corners2color mean

splot '../../Data/Structure/kartesianDeformations_0' u 1:2:4 notitle with pm3d,\
      '../../Data/Structure/kartesianDeformations_1' u 1:2:4 notitle with pm3d,\
      '../../Data/Structure/kartesianDeformations_2' u 1:2:4 notitle with pm3d,\
      '../../Data/Structure/kartesianDeformations_3' u 1:2:4 notitle with pm3d