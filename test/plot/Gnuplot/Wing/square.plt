set terminal png background rgb 'black' size 1000, 1000
set output '../../png/Wing/square.png'

set palette defined (0 "blue", 0.5 "green", 0.8 "yellow", 1 "red")
set cbrange [0:0.1]
set xlabel 'x' tc rgb 'gray'
set ylabel 'y' tc rgb 'gray'
set zlabel 'dcp' tc rgb 'gray'
set key tc rgb 'gray'
set border lc rgb 'gray'
set pm3d map
# set pm3d interpolate 10,10 corners2color mean

set multiplot layout 2,2 rowsfirst
set title "Rotated by 0°" tc rgb 'gray'
splot '../../Data/Wing/square0' notitle with pm3d
set title "Rotated by 90°" tc rgb 'gray'
splot '../../Data/Wing/square1' notitle with pm3d
set title "Rotated by 180°" tc rgb 'gray'
splot '../../Data/Wing/square2' notitle with pm3d
set title "Rotated by 270°" tc rgb 'gray'
splot '../../Data/Wing/square3' notitle with pm3d
unset multiplot