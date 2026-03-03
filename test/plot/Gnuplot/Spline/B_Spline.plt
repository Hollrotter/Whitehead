set terminal png background rgb 'black' size 2000, 500

set xlabel 'x' tc rgb 'gray'
set ylabel 'y' tc rgb 'gray'
set zlabel 'z' tc rgb 'gray'
set border lc rgb 'gray'
set grid
set key tc rgb 'gray'
set key opaque
set ytics 0,0.2,1
set output '../../png/Spline/B_Spline.png'

plot '../../Data/Spline/B_Spline' u 1:2   with lines notitle lw 3 lc rgb 'white',\
     '../../Data/Spline/B_Spline' u 1:3   with lines notitle lw 3 lc rgb 'white',\
     '../../Data/Spline/B_Spline' u 1:4   with lines notitle lw 3 lc rgb 'white',\
     '../../Data/Spline/B_Spline' u 1:5   with lines notitle lw 3 lc rgb 'white',\
     '../../Data/Spline/B_Spline' u 1:6   with lines notitle lw 3 lc rgb 'white',\
     '../../Data/Spline/B_Spline' u 1:7   with lines notitle lw 3 lc rgb 'white',\
     '../../Data/Spline/B_Spline' u 1:8   with lines notitle lw 3 lc rgb 'white',\
     '../../Data/Spline/B_Spline' u 1:9   with lines notitle lw 3 lc rgb 'white',\
     '../../Data/Spline/B_Spline' u 1:10  with lines notitle lw 3 lc rgb 'white',\
     '../../Data/Spline/B_Spline' u 1:11  with lines notitle lw 3 lc rgb 'white',\
     '../../Data/Spline/B_Spline' u 1:12  with lines notitle lw 3 lc rgb 'white',\
     '../../Data/Spline/B_Spline' u 1:13  with lines notitle lw 3 lc rgb 'white'