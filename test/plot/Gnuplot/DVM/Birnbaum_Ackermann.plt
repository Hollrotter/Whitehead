set terminal png size 500, 250
# set terminal png background rgb 'black' size 2000, 1000
# set border lc rgb 'gray'
# set key tc rgb 'gray'
set grid
set xlabel 'x/m'
set ylabel 'gamma/m'
set yrange [0:2]
set output '../../png/DVM/Birnbaum_Ackermann.png'

plot '../../Data/DVM/dvm1'      notitle with lines lw 3 lc rgb 'red',\
     '../../Data/DVM/dcp' u 1:2 notitle with lines lw 3 lc rgb 'blue'