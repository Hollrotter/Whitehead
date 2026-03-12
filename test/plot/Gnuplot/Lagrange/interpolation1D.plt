set terminal png size 500, 200

set xlabel 'x'
set ylabel 'z' offset 2
set yrange [-4:10]
set xtics -1,0.5,1
set ytics -4,2,10
set grid
set output '../../png/Lagrange/interpolation1D.png'

plot '../../Data/Lagrange/interpolation1D_1' notitle with points pt 7 ps 1,\
     '../../Data/Lagrange/interpolation1D_2' notitle with lines lw 1 lc rgb 'red'