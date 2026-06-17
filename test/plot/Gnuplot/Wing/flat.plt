set terminal png background rgb 'black' size 2000, 1000
set output '../../png/Wing/flat.png'

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
set cbrange [0:]
set ticslevel 0.0
set xlabel 'y' tc rgb 'gray'
set ylabel 'x' tc rgb 'gray'
set zlabel 'dcp' tc rgb 'gray'
set key tc rgb 'gray'
set border lc rgb 'gray'
set pm3d map
# set pm3d interpolate 10,10 corners2color mean

splot '../../Data/Wing/flat' u 2:1:5 notitle with pm3d