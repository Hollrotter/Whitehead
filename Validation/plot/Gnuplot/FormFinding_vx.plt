set terminal png size 600, 450

set output '../png/FormFinding_vx.png'

# set xrange [-.5:.5]
# set yrange [-.5:.5]
# set zrange [0:]
# set palette defined (0 "black", 0.01 "#000080", 0.25 "#0080ff", 1 "white")
# set palette defined (0 "blue", 0.5 "green", 0.8 "yellow", 1 "red")
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
set zlabel 'v_x/m'

set ticslevel 0.0
set xtics -30,10,30
set ytics -10, 5,10

set zlabel offset 1.0,0.0
set xtics offset 0,-0.5

set pm3d interpolate 10,10 corners2color mean

splot '../Data/FormFinding_vx_0' u 1:2:3 notitle with pm3d,\
      '../Data/FormFinding_vx_1' u 1:2:3 notitle with pm3d