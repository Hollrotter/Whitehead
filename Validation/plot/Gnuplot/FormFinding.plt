set terminal png size 1000, 1000

set output '../png/FormFinding.png'

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

set pm3d interpolate 10,10 corners2color mean

splot '../Data/FormFinding' notitle with pm3d

# splot '../Data/FormFinding_0' notitle with pm3d,\
#       '../Data/FormFinding_1' notitle with pm3d,\
#       '../Data/FormFinding_2' notitle with pm3d,\
#       '../Data/FormFinding_3' notitle with pm3d,\
#       '../Data/FormFinding_4' notitle with pm3d,\
#       '../Data/FormFinding_5' notitle with pm3d,\
#       '../Data/FormFinding_6' notitle with pm3d,\
#       '../Data/FormFinding_7' notitle with pm3d,\
#       '../Data/FormFinding_8' notitle with pm3d