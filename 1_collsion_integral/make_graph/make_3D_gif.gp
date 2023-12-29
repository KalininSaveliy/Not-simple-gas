count = 400
step = 1

set output "graph/result.gif"
set title "Relaxation"

set xlabel 'Vx'
set ylabel 'Vy'
set zlabel 'f'

set xrange [-5:5]
set yrange [-5:5]
# set zrange [0:0.07]
set zrange [0 : 0.035]


set grid
set pm3d
# set cbrange [0:0.07]
set cbrange [0 : 0.035]
set ticslevel 0  # 0 оси Z на плоскости XY
set view 30,30
# set view map

set term gif animate delay 1 size 1000, 800
# set term gif animate delay 1 size 1500, 800

do for [n = 0 : count : step] {
    if (n == 0) {
        filename = sprintf("data/init.txt")
    } else {
        filename = sprintf("data/f_%d.txt", n)
    }
    splot filename w pm3d title "Velocity"
}
