set terminal png background rgb 'black' size 1200, 1000

set output '../../png/Aerodynamics/panel.png'

# set palette defined (0 "black", 0.01 "#000080", 0.25 "#0080ff", 1 "white")
set palette defined (0  0.0 0.0 1.0, \
                     1  0.0 0.698 1.0, \
                     2  0.0 1.0 1.0, \
                     3  0.0 1.0 0.698, \
                     4  0.0 1.0 0.0, \
                     5  0.698 1.0 0.0, \
                     6  1.0 1.0 0.0, \
                     7  1.0 0.698 0.0, \
                     8  1.0 0.0 0.0 )
# set palette maxcolor 100
set zrange [0:]
set cbrange [0:]
set ticslevel 0.0
set xlabel 'y' tc rgb 'gray'
set ylabel 'x' tc rgb 'gray'
set zlabel 'dcp' tc rgb 'gray'
set key tc rgb 'gray'
set border lc rgb 'gray'

set view map
# set pm3d depthorder base
# set pm3d interpolate 10,10 corners2color mean

splot '../../Data/Aerodynamics/panel_0'  u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_1'  u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_2'  u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_3'  u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_4'  u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_5'  u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_6'  u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_7'  u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_8'  u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_9'  u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_10' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_11' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_12' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_13' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_14' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_15' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_16' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_17' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_18' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_19' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_20' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_21' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_22' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_23' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_24' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_25' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_26' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_27' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_28' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_29' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_30' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_31' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_32' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_33' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_34' u 2:1:4 notitle with pm3d,\
      '../../Data/Aerodynamics/panel_35' u 2:1:4 notitle with pm3d