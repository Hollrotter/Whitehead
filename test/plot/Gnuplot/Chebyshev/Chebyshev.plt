set terminal png background rgb 'black' size 1500, 750

set border lc rgb 'gray'
set key tc rgb 'gray'
set xlabel 'x' tc rgb 'gray'
set ylabel 'z' tc rgb 'gray'
set grid
set output '../../png/Chebyshev.png'

plot '../../Data/Chebyshev' u 1:2 notitle with lines lw 4 lc rgb 'red',\
     '../../Data/Chebyshev' u 1:3 notitle with lines lw 4 lc rgb 'blue',\
     '../../Data/Chebyshev' u 1:4 notitle with lines lw 4 lc rgb 'green'