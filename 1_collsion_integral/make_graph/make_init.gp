count = 10000
step = 50

set output "graph/init.jpg"
set title "Relaxation"

set xlabel 'Vx'
set ylabel 'Vy'
set zlabel 'f'

set xrange [-5:5]
set yrange [-5:5]
set zrange [0:0.035]

set grid
set pm3d
set cbrange [0:0.035]
set ticslevel 0  # 0 оси Z на плоскости XY
set view 30,210
# set view map

# set size 1500, 800

filename = "data/init.txt"

splot filename w pm3d title "Velocity"