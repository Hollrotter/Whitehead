set terminal png background rgb 'black' size 2000, 1000
set output '../../png/VLM/vlm.png'

set palette defined (0 "blue", 0.5 "green", 0.8 "yellow", 1 "red")
set zrange [0:2.5]
set cbrange [0:2.5]
set xlabel 'y' tc rgb 'gray'
set ylabel 'x' tc rgb 'gray'
set zlabel 'dcp' tc rgb 'gray'
set key tc rgb 'gray'
set border lc rgb 'gray'
# set pm3d map

splot '../../Data/VLM/vlm' u 2:1:3 notitle with pm3d