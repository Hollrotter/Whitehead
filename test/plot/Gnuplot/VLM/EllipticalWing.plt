set terminal png enhanced size 750, 400
# set terminal png background rgb 'black' size 2000, 1000
set output '../../png/VLM/EllipticalWing.png'

# set palette defined (0 "blue", 0.5 "green", 0.8 "yellow", 1 "red")
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
# set cbrange [-0.8:0.6]
# set zrange [-1:]
set xlabel 'y/m'# tc rgb 'gray'
set ylabel 'x/m'# tc rgb 'gray'
# set key tc rgb 'gray'
# set border lc rgb 'gray'
set view map
set pm3d interpolate 10,10 corners2color mean

splot '../../Data/VLM/EllipticalWing' u 2:1:3 notitle with pm3d