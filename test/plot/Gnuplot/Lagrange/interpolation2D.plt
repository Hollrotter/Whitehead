set terminal png background rgb 'black' size 1000, 1000

set output '../../png/Lagrange/interpolation2D.png'
# set cbrange [-3:3]
set palette defined (0 "black", 0.337 "#000080", 0.675 "#0080ff", 1 "white")
set palette maxcolor 100
set xlabel 'x' tc rgb 'gray'
set ylabel 'y' tc rgb 'gray'
set key tc rgb 'gray'
set border lc rgb 'gray'

set pm3d map
set pm3d interpolate 10,10 corners2color mean

splot '../../Data/Lagrange/interpolation2D' notitle