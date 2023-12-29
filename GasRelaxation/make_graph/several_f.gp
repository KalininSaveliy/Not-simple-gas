# set terminal png size 800, 600
set output ".\\graph\\several_f.png"

set title "Relaxation"

set xlabel "t (во временах свободного пробега)"
set ylabel "M_{xx}"

# set logscale y
set yrange [2:5]

# set format y '%.1e'


set xtics 5
set mxtics 5
set ytics 0.5
set mytics 5
set style line 100 lt 1 lc rgb "dark-gray" lw 1  # for major
set style line 101 lt 0.5 lc rgb "gray" lw 1  # for minor
set grid mytics ytics ls 100, ls 101
set grid mxtics xtics ls 100, ls 101

name = ".\\data\\several_f.txt"
stats name
n_curv = STATS_columns - 1
array Names[n_curv]
Names[1] = sprintf("1: N_v = %d, n_{col} = %d", 10, 1000)
Names[2] = sprintf("2: N_v = %d, n_{col} = %d", 10, 5000)
Names[3] = sprintf("3: N_v = %d, n_{col} = %d", 10, 10000)
Names[4] = sprintf("4: N_v = %d, n_{col} = %d", 10, 25000)
Names[5] = sprintf("5: N_v = %d, n_{col} = %d", 30, 250000)

plot\
     for [i = 1:n_curv] name u 1:(column(i+1)) title Names[i] w l
pause -1
#  name u 1:2 title "M_{1 xx} h = 2.50" w l lc rgb "red",\
#      name u 1:3 title "M_{2 xx} h = 1.67" w l lc rgb "green",\
#      name u 1:4 title "M_{3 xx} h = 1.25" w l lc rgb "violet",\
#      name u 1:5 title "M_{4 xx} h = 1.00" w l lc rgb "blue",\
#      name u 1:6 title "M_{5 xx} h = 0.33" w l lc rgb "orange"
