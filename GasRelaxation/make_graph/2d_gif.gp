set output "graph/sample_2D.gif"
set title "Relaxation (n_v = 30, n_{col} = 250000, \\tau = 0.02)"

set xlabel'v, в скоростях sqrt(k * T_0 / m)'
set ylabel 'f(v)'

set xrange [-5:5]
set yrange [0 : 0.3]

set xtics 1
set mytics 10
set mxtics 10

set style line 100 lt 1 lc rgb "dark-gray" lw 1  # for major
set style line 101 lt 0.5 lc rgb "gray" lw 1  # for minor
set grid mytics ytics ls 100, ls 101
set grid mxtics xtics ls 100, ls 101


set term gif animate delay 1 size 1000, 800

count = 1000
step = 1

do for [n = 0: count : step] {
    filename = sprintf("data/sample_f/f%d.txt", n)
    plot filename w lp title sprintf("n_{iter} = %d / %d", n, count)
}
