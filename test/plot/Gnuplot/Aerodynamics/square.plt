set terminal png background rgb 'black' size 2000, 1000
set output '../../png/Aerodynamics/square.png'

set palette defined (0 "blue", 0.5 "green", 0.8 "yellow", 1 "red")
# set cbrange [0:0.2]
set xlabel 'x' tc rgb 'gray'
set ylabel 'y' tc rgb 'gray'
set zlabel 'dcp' tc rgb 'gray'
set key tc rgb 'gray'
set border lc rgb 'gray'

set pm3d map
# set pm3d interpolate 10,10 corners2color mean

z = 4
set multiplot layout 4,4 rowsfirst
set title "w1 rotated by 0°, w2 rotated by 0°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square0_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square0_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 0°, w2 rotated by 90°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square1_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square1_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 0°, w2 rotated by 180°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square2_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square2_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 0°, w2 rotated by 270°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square3_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square3_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 90°, w2 rotated by 0°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square4_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square4_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 90°, w2 rotated by 90°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square5_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square5_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 90°, w2 rotated by 180°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square6_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square6_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 90°, w2 rotated by 270°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square7_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square7_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 180°, w2 rotated by 0°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square8_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square8_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 180°, w2 rotated by 90°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square9_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square9_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 180°, w2 rotated by 180°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square10_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square10_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 180°, w2 rotated by 270°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square11_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square11_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 270°, w2 rotated by 0°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square12_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square12_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 270°, w2 rotated by 90°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square13_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square13_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 270°, w2 rotated by 180°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square14_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square14_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 270°, w2 rotated by 270°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square15_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square15_1' u 1:2:z notitle with pm3d
unset multiplot