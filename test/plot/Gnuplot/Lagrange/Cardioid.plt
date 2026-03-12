set terminal png background rgb 'black' size 1000, 1000

set border lc rgb 'gray'
set key tc rgb 'gray'
set xlabel 'x' tc rgb 'gray'
set ylabel 'y' tc rgb 'gray'
# set xtics 0,0.5,2
# set ytics 0,0.25,2
set grid
set output '../../png/Lagrange/Cardioid.png'

plot '../../Data/Lagrange/Cardioid' u 1:2 notitle with lines lw 1 lc rgb 'red',\
     '../../Data/Lagrange/Cardioid' u 3:4 notitle pointsize 2