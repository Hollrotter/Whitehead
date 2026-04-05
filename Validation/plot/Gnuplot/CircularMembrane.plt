set terminal png transparent size 610, 380
# set terminal svg size 750, 500

# set format y "%.1f"
set xlabel 'x/a'
set ylabel 'z/a' offset 0,-7
set yrange [0:0.6]
set xtics 0,0.1,1
set ytics 0,0.2,0.6
set format x ""
set format y ""

set grid

set output '../png/CircularMembrane.png'

plot '../Data/CircularMembrane_0.txt' notitle with points pt 7 ps 2,\
     '../Data/CircularMembrane_1.txt' notitle with points pt 7 ps 2,\
     '../Data/CircularMembrane_2.txt' notitle with points pt 7 ps 2,\
     '../Data/CircularMembrane_3.txt' notitle with points pt 7 ps 2,\
     '../Data/CircularMembrane_4.txt' notitle with points pt 7 ps 2,\
     '../Data/CircularMembrane_5.txt' notitle with points pt 7 ps 2,\
     '../Data/CircularMembrane_6.txt' notitle with points pt 7 ps 2,\
     '../Data/CircularMembrane_7.txt' notitle with points pt 7 ps 2