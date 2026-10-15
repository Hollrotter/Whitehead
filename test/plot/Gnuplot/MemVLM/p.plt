set terminal png size 500, 450

set output '../../png/MemVLM/p.png'
set cbrange [0:0.9]
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

set xtics -1,0.25
set ytics  0,0.25
set cbtics 0,0.1
set lmargin 6
set xlabel 'x/m'
set ylabel 'y/m'

set view map
set pm3d interpolate 10,10 corners2color mean

splot '../../Data/MemVLM/p' notitle with pm3d