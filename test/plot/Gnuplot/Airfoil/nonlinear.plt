set terminal png background rgb 'black' size 1000, 1000

set xrange [0:2]
# set yrange [0:5]
set border lc rgb 'gray'
set key tc rgb 'gray'
set xlabel 'x'   tc rgb 'gray'
set ylabel 'dcp' tc rgb 'gray'
set grid
set output '../../png/Airfoil/nonlinear.png'

plot '../../Data/Airfoil/nonlinear' notitle with lines lw 4 lc rgb 'red'