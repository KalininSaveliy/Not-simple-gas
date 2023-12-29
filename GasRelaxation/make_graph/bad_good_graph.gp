# set terminal png
set output ".\\graph\\bad_good.png"

set title "Relaxation"

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

name = ".\\data\\bad and good.txt"
plot name u 1:2 title "M_{bad xx}" w l lc rgb "red",\
     name u 1:3 title "M_{bad yy}" w l lc rgb "green",\
     name u 1:4 title "M_{bad zz}" w l lc rgb "violet",\
     name u 1:5 title "M_{good xx}" w l lc rgb "blue",\
     name u 1:6 title "M_{good yy}" w l lc rgb "orange",\
     name u 1:7 title "M_{good zz}" w l lc rgb "dark-plum"
