set terminal png background rgb 'black' size 1500, 750

set border lc rgb 'gray'
set key tc rgb 'gray'
set xlabel 'x' tc rgb 'gray'
set ylabel 'z' tc rgb 'gray'
set grid
set output '../../png/Lagrange/derivativeMatrix.png'

plot '../../Data/Lagrange/derivativeMatrix' u 1:2 notitle with lines lw 4 lc rgb 'red',\
     '../../Data/Lagrange/derivativeMatrix' u 1:3 notitle with lines lw 4 lc rgb 'blue'