set terminal png size 500, 400

set xlabel 'x'
set ylabel 'z' offset 2
set xtics -1,0.5,1
set grid
set output '../../png/Lagrange/interpolation1D.png'

set multiplot layout 2,1 rowsfirst
set yrange [-4:10]
set ytics -4,2,10
set lmargin 6.5
plot '../../Data/Lagrange/interpolation1D_1' notitle with points pt 7 ps 1,\
     '../../Data/Lagrange/interpolation1D_2' notitle with lines lw 1 lc rgb 'red'
set yrange [0:0.12]
set ytics 0,0.04,0.12
set lmargin 6.5
plot '../../Data/Lagrange/interpolation1D_3' notitle with lines lw 1 lc rgb 'blue'
unset multiplot