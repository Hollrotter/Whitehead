set terminal png background rgb 'black' size 1000, 1000

set border lc rgb 'gray'
set key tc rgb 'gray'
set xlabel 'x'   tc rgb 'gray'
set ylabel 'dcp' tc rgb 'gray'
set grid
set output '../../png/Jacobi/JacobiB.png'

plot '../../Data/Jacobi/JacobiB' u 1:2 notitle with lines lw 2 lc rgb 'red',\
     '../../Data/Jacobi/JacobiB' u 1:3 notitle with lines lw 2 lc rgb 'blue',\
     '../../Data/Jacobi/JacobiB' u 1:4 notitle with lines lw 2 lc rgb 'green',\
     '../../Data/Jacobi/JacobiB' u 1:5 notitle with lines lw 2 lc rgb 'yellow',\
     '../../Data/Jacobi/JacobiB' u 1:6 notitle with lines lw 2 lc rgb 'brown',\
     '../../Data/Jacobi/JacobiB' u 1:7 notitle with lines lw 2 lc rgb 'white',\
     '../../Data/Jacobi/JacobiB' u 1:8 notitle with lines lw 2 lc rgb 'violett'