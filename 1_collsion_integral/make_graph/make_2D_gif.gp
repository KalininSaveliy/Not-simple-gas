set output "graph/result_2D.gif"
set title "Relaxation (n_v = 30; v_y, v_z = const)"

set xlabel'v, в скоростях sqrt(k * T_0 / m)'
set ylabel 'f(v)'

set xrange [-5:5]
set yrange [0 : 0.035]

set mytics 10
set mxtics 10
set grid mxtics mytics

set term gif animate delay 1 size 1000, 800

count = 300
step = 1

do for [n = 0: count : step] {
    if (n == 0) {
        filename = sprintf("data/init.txt")
    } else {
        filename = sprintf("data/f_%d.txt", n)
    }
    plot filename w lp title sprintf("n_{iter} = %d / %d", n, count)
}
