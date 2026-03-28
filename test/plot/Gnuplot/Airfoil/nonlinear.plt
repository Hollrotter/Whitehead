set terminal png background rgb 'black' size 1000, 1000

set xrange [0:1]
set yrange [:2]
set border lc rgb 'gray'
set key tc rgb 'gray'
set xlabel 'x'   tc rgb 'gray'
set ylabel 'dcp' tc rgb 'gray'
set xtics 0,0.2,1
set ytics 0,0.25,2
set grid
set output '../../png/Airfoil/nonlinear.png'

plot '../../Data/Airfoil/nonlinear' notitle with lines lw 4 lc rgb 'red'