set terminal png background rgb 'black' size 2000, 1000

set output '../../png/Aerodynamics/vMesh.png'

set xlabel 'y/m' tc rgb 'gray'
set ylabel 'x/m' tc rgb 'gray'
set key tc rgb 'gray'
set border lc rgb 'gray'

set view map

splot '../../Data/Aerodynamics/v_0' u 2:1:4 notitle with lines,\
      '../../Data/Aerodynamics/v_1' u 2:1:4 notitle with lines,\
      '../../Data/Aerodynamics/v_2' u 2:1:4 notitle with lines