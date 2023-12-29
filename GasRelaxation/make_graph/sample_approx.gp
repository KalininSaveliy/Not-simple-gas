set terminal png size 800, 600

filename = ".\\data\\sample_m.txt"
set output ".\\graph\\sample_approx.png"
set title "Approx"

set xlabel "t"
set ylabel sprintf("log(M_2 + m)")

# set xrange[0:20]
# set logscale y
# set yrange [1:5]
# set format y '%.1e'

set xtics 5
set mxtics 5
set ytics 2
set mytics 10
set style line 100 lt 1 lc rgb "dark-gray" lw 1  # for major
set style line 101 lt 0.5 lc rgb "gray" lw 1  # for minor
set grid mytics ytics ls 100, ls 101
set grid mxtics xtics ls 100, ls 101


n_curv = 3
array Shifts[n_curv]
stats filename u 2 name "X"
stats filename u 3 name "Y"
stats filename u 4 name "Z"
Shifts[1] = -X_min
Shifts[2] = -Y_max
Shifts[3] = -Z_max



# fitting with polynom
array K[n_curv]
array B[n_curv]
x_l = 1
y_l = -8
set label 1 "y = k * x + b" at x_l, y_l
f(x) = k * x + b
fit [0.0:15.0] f(x) filename u 1:(log(column(2) + Shifts[1])) via k, b
set label 2 sprintf("%d: k = %f, b = %f", 1, k, b) at x_l, y_l - 0.6

do for[i = 2:n_curv] {
    fit [0.0:15.0] f(x) filename u 1:(log(-(column(i+1) + Shifts[i]))) via k, b
    set label i+1 sprintf("%d: k = %f, b = %f", i, k, b) at x_l, y_l - 0.6 * i
    K[i] = k
    B[i] = b
}
# end of fitting

array Names[n_curv]
Names[1] = "log(M_{xx} + m_x)"
Names[2] = "log(M_{yy} + m_y)"
Names[3] = "log(M_{zz} + m_z)"

plot filename u 1:(log(column(2) + Shifts[1])) title Names[1] w l,\
     filename u 1:(log(-(column(3) + Shifts[2]))) title Names[2] w l,\
     filename u 1:(log(-(column(4) + Shifts[3]))) title Names[3] w l

# plot\
#     for [i = 2:n_curv] filename u 1:(log(column(i+1) + Shifts[i])) title Names[i] w l

pause -1  # процесс не завершен, поэтому работает зуууум и прочее
