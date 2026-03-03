set terminal png background rgb 'black' size 2000, 500

set xlabel 'x' tc rgb 'gray'
set ylabel 'y' tc rgb 'gray'
set border lc rgb 'gray'
set grid
set key tc rgb 'gray'
set key opaque
set ytics 0,0.2,1
set output '../../png/Spline/Splinefit.png'

plot '../../Data/Spline/Splinefit' u 1:2 with lines notitle linewidth 4 linecolor rgb "red",\
     '../../Data/Spline/Splinefit' u 1:3 notitle