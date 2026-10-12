set terminal png background rgb 'black' size 500, 250

set xrange [0:1]
set yrange [0:2]
set border lc rgb 'gray'
set key tc rgb 'gray'
set xlabel 'x'   tc rgb 'gray'
set ylabel 'dcp' tc rgb 'gray'
set grid
set output '../../png/MemDVM/p.png'

plot '../../Data/MemDVM/p' notitle with lines lw 4 lc rgb 'red'