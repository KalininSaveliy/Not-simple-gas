# set terminal png
set output "graph/time_1000003.png"

set title "Relaxation time"

set xlabel "t (во временах свободного пробега)"
set ylabel "M_{2}"

# set logscale y
set yrange [1:5]

# set format y '%.1e'


set xtics 5
set mxtics 5
set ytics 0.5
# set mytics 1
set grid xtics ytics mxtics mytics

name = "data/relax_time_250007.txt"
plot name u 1:2 title "M_{xx}" w l lc rgb "red",\
     name u 1:3 title "M_{yy}" w l lc rgb "green",\
     name u 1:4 title "M_{zz}" w l lc rgb "blue"
