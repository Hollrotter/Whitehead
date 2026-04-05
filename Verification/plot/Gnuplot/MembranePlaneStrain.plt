set terminal png background rgb 'black' size 1200, 1000

# set xrange [0:1]
# set xtics 0,0.2,1
# set yrange [0:8]
# set ytics 0,1,8
set format y "%.2e"
set border lc rgb 'gray'
set key tc rgb 'gray'
set logscale y 10
# set xlabel 'N' tc rgb 'gray'
# set ylabel 'U/Z' tc rgb 'gray'
set grid
set output '../png/MembranePlaneStrainMMS.png'

plot '../Data/MembranePlaneStrainMMS' u 1:3 notitle with lines lw 1 lc rgb 'red'