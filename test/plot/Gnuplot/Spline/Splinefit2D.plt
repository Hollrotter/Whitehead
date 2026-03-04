set terminal png background rgb 'black' size 2000, 500

set xlabel 'x' tc rgb 'gray'
set ylabel 'y' tc rgb 'gray'
set zlabel 'z' tc rgb 'gray'
set border lc rgb 'gray'
set grid
set key tc rgb 'gray'
set key opaque
set output '../../png/Spline/Splinefit2D.png'

splot '../../Data/Spline/Splinefit2D' u 1:2:5 notitle with points pt 7, \
      '../../Data/Spline/z'           u 1:2:5 notitle with lines