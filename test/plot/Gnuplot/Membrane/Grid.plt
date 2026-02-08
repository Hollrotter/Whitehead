set terminal png size 1000, 1000

set output '../../png/Membrane/xy.png'

set xlabel 'x'
set ylabel 'y'

# unset border
# unset xtics
# unset ytics

set view map

splot '../../Data/Membrane/xy' notitle with lines lw 2 lc rgb 'black'