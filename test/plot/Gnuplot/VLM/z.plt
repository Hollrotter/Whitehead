set terminal png background rgb 'black' size 2000, 1000
set output '../../png/VLM/z.png'

set palette defined (0 "blue", 0.5 "green", 0.8 "yellow", 1 "red")
# set cbrange [-0.5:1.0]
set xlabel 'y' tc rgb 'gray'
set ylabel 'x' tc rgb 'gray'
set zlabel 'dcp' tc rgb 'gray'
set key tc rgb 'gray'
set border lc rgb 'gray'
set pm3d map

splot '../../Data/VLM/z' u 2:1:4 notitle