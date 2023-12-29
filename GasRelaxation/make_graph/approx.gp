# set terminal png size 800, 600

filename = ".\\data\\sample_m.txt"
set output ".\\graph\\approx_time.png"
set title "Approx"

set xlabel "t"
set ylabel sprintf("log(M_{xx} + shift)")

# set xrange[0:20]
# set logscale y
# set yrange [1:5]
# set format y '%.1e'

set xtics 5
set mxtics 5
set ytics 0.5
# set mytics 1
set grid xtics ytics mxtics mytics


stats filename
n_curv = STATS_columns - 1
array Shifts[n_curv]
do for[i = 1:n_curv] {
    stats filename u i+1 name "A"
    Shifts[i] = -A_min
}

# fitting with polynom
array K[n_curv]
array B[n_curv]
x_l = 1
y_l = -8
set label 1 "y = k * x + b" at x_l, y_l
f(x) = k * x + b
do for[i = 1:n_curv] {
    fit [0.0:15.0] f(x) filename u 1:(log(column(i+1) + Shifts[i])) via k, b
    set label i+1 sprintf("%d: k = %f, b = %f, tau = %f", i, k, b, -1/k) at x_l, y_l - 0.5 * i
    K[i] = k
    B[i] = b
}
# end of fitting

array Names[n_curv]
Names[1] = sprintf("1: N_v = %d, n_{col} = %d", 10, 1000)
Names[2] = sprintf("2: N_v = %d, n_{col} = %d", 10, 5000)
Names[3] = sprintf("3: N_v = %d, n_{col} = %d", 10, 10000)
Names[4] = sprintf("4: N_v = %d, n_{col} = %d", 10, 25000)
Names[5] = sprintf("5: N_v = %d, n_{col} = %d", 30, 250000)

plot\
    for [i = 1:n_curv] filename u 1:(log(column(i+1) + Shifts[i])) title Names[i] w l

pause -1  # процесс не завершен, поэтому работает зуууум и прочее
