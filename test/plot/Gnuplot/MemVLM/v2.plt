set terminal png background rgb 'black' size 1000, 1000

set output '../../png/MemVLM/v2.png'
# set cbrange [0:0.2]
set palette defined (0 "black", 0.01 "#000080", 0.25 "#0080ff", 1 "white")
set palette maxcolor 100

set xlabel 'x/m' tc rgb 'gray'
set ylabel 'y/m' tc rgb 'gray'
set key tc rgb 'gray'
set border lc rgb 'gray'

set view map
set pm3d interpolate 10,10 corners2color mean

splot '../../Data/MemVLM/v2' notitle with pm3d