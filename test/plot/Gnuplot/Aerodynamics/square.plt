set terminal png background rgb 'black' size 2000, 1000
set output '../../png/Aerodynamics/square.png'

set palette defined (0 "blue", 0.5 "green", 0.8 "yellow", 1 "red")
set cbrange [0:0.18]
set xlabel 'x' tc rgb 'gray'
set ylabel 'y' tc rgb 'gray'
set key tc rgb 'gray'
set border lc rgb 'gray'

set pm3d map
# set pm3d interpolate 10,10 corners2color mean

z = 4
set multiplot layout 4,4 rowsfirst
set title "w1 rotated by 0°, w2 rotated by 0°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square00_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square00_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 0°, w2 rotated by 90°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square01_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square01_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 0°, w2 rotated by 180°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square02_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square02_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 0°, w2 rotated by 270°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square03_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square03_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 90°, w2 rotated by 0°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square10_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square10_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 90°, w2 rotated by 90°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square11_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square11_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 90°, w2 rotated by 180°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square12_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square12_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 90°, w2 rotated by 270°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square13_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square13_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 180°, w2 rotated by 0°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square20_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square20_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 180°, w2 rotated by 90°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square21_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square21_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 180°, w2 rotated by 180°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square22_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square22_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 180°, w2 rotated by 270°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square23_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square23_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 270°, w2 rotated by 0°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square30_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square30_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 270°, w2 rotated by 90°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square31_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square31_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 270°, w2 rotated by 180°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square32_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square32_1' u 1:2:z notitle with pm3d
set title "w1 rotated by 270°, w2 rotated by 270°" tc rgb 'gray'
splot '../../Data/Aerodynamics/square33_0' u 1:2:z notitle with pm3d,\
      '../../Data/Aerodynamics/square33_1' u 1:2:z notitle with pm3d
unset multiplot