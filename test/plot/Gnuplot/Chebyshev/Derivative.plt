set terminal png background rgb 'black' size 750, 750

set border lc rgb 'gray'
set key tc rgb 'gray'
set xlabel 'x' tc rgb 'gray'
set ylabel 'z' tc rgb 'gray'
set grid
set output '../../png/Chebyshev/Derivative.png'

plot '../../Data/Chebyshev/Derivative' u 1:2 notitle with lines lw 4 lc rgb 'red',\
     '../../Data/Chebyshev/Derivative' u 1:3 notitle with lines lw 4 lc rgb 'blue',\
     '../../Data/Chebyshev/Derivative' u 1:4 notitle with lines lw 4 lc rgb 'green'