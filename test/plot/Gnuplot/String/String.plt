set terminal png size 500, 250
# set terminal png background rgb 'black' size 2000, 1000
# set border lc rgb 'gray'
# set key tc rgb 'gray'
set grid
set xlabel 'x/m'
set ylabel 'z/m'
set output '../../png/String/Linear.png'

plot '../../Data/String/Linear' notitle with lines lw 3 lc rgb 'red'