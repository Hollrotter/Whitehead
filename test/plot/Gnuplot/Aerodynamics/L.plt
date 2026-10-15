set terminal wxt background rgb 'black' size 1000, 1000

# set output '../../png/Aerodynamics/L.png'
# set cbrange [0:]
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

set xlabel 'y/m' tc rgb 'gray'
set ylabel 'x/m' tc rgb 'gray'
set zlabel 'dcp' tc rgb 'gray'
set key tc rgb 'gray'
set border lc rgb 'gray'

# set view map
set pm3d depthorder base
# set pm3d interpolate 10,10 corners2color mean

splot '../../Data/Aerodynamics/L_0' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/L_1' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/L_2' u 2:1:4 notitle with pm3d