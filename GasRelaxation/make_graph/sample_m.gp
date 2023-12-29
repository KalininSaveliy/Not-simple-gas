set terminal png size 800, 600
set output ".\\graph\\sample_m.png"

set title "Relaxation"

set xlabel "t"
set ylabel "M_2"

# set logscale y
set yrange [1:5]

# set format y '%.1e'
set xtics 5
set mxtics 5
set ytics 0.5
set mytics 5
set style line 100 lt 1 lc rgb "dark-gray" lw 1  # for major
set style line 101 lt 0.5 lc rgb "gray" lw 1  # for minor
set grid mytics ytics ls 100, ls 101
set grid mxtics xtics ls 100, ls 101

name = ".\\data\\sample_m.txt"

plot name u 1:2 title "M_{xx}" w l lc rgb "red",\
     name u 1:3 title "M_{yy}" w l lc rgb "green",\
     name u 1:4 title "M_{zz}" w l lc rgb "blue"
pause -1
