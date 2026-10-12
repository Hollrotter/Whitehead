set terminal png background rgb 'black' size 500, 250

set xrange [0:1]
set yrange [0:0.025]
set border lc rgb 'gray'
set key tc rgb 'gray'
set xlabel 'x' tc rgb 'gray'
set ylabel 'z' tc rgb 'gray'
set grid
set output '../../png/MemDVM/z.png'

plot '../../Data/MemDVM/z' notitle with lines lw 3 lc rgb 'red'