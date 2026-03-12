set terminal png background rgb 'black' size 2000, 1000

set output '../../png/Lagrange/interpolationCurvilinear.png'
# set cbrange [-3:3]
set palette defined (0 "black", 0.337 "#000080", 0.675 "#0080ff", 1 "white")
set palette maxcolor 100
set xlabel 'y' tc rgb 'gray'
set ylabel 'x' tc rgb 'gray'
set key tc rgb 'gray'
set border lc rgb 'gray'

set pm3d map
set pm3d interpolate 10,10 corners2color mean

splot '../../Data/Lagrange/interpolationCurvilinear' u 2:1:3 notitle