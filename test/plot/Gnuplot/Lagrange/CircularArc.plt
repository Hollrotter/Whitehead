set terminal png background rgb 'black' size 1500, 750

set border lc rgb 'gray'
set key tc rgb 'gray'
set xlabel 'x' tc rgb 'gray'
set ylabel 'y' tc rgb 'gray'
set xrange[-1:1]
set yrange[-1:1]
set grid
set output '../../png/Lagrange/CircularArc.png'

set multiplot layout 2,4 rowsfirst
set title "0°-90°" tc rgb 'gray'
plot '../../Data/Lagrange/CircularArc' u  1:2  notitle with lines lw 2
set title "45°-135°" tc rgb 'gray'
plot '../../Data/Lagrange/CircularArc' u  3:4  notitle with lines lw 2
set title "90°-180°" tc rgb 'gray'
plot '../../Data/Lagrange/CircularArc' u  5:6  notitle with lines lw 2
set title "135°-225°" tc rgb 'gray'
plot '../../Data/Lagrange/CircularArc' u  7:8  notitle with lines lw 2
set title "180°-270°" tc rgb 'gray'
plot '../../Data/Lagrange/CircularArc' u  9:10 notitle with lines lw 2
set title "225°-315°" tc rgb 'gray'
plot '../../Data/Lagrange/CircularArc' u 11:12 notitle with lines lw 2
set title "270°-360°" tc rgb 'gray'
plot '../../Data/Lagrange/CircularArc' u 13:14 notitle with lines lw 2
set title "315°-45°" tc rgb 'gray'
plot '../../Data/Lagrange/CircularArc' u 15:16 notitle with lines lw 2
unset multiplot