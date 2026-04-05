set terminal png size 1500, 500

set output '../png/Kirsch_Kt.png'
set xrange [0:0.75]
set yrange [0:0.25]
set cbrange [-0.25:2.25]
set cbtics -0.25,0.25,2.25
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

set view map
set pm3d interpolate 10,10 corners2color mean

splot '../Data/Kirsch_s1_0' u 1:2:($3/33333.333333333) notitle with pm3d,\
      '../Data/Kirsch_s1_1' u 1:2:($3/33333.333333333) notitle with pm3d,\
      '../Data/Kirsch_s1_2' u 1:2:($3/33333.333333333) notitle with pm3d,\
      '../Data/Kirsch_s1_3' u 1:2:($3/33333.333333333) notitle with pm3d,\
      '../Data/Kirsch_s1_4' u 1:2:($3/33333.333333333) notitle with pm3d