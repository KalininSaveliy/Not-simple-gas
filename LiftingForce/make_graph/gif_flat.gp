count = 100
step = 10
max_z = 1.2

# set palette defined (0 "white", max_z "red")
set cbrange [0:max_z]

set xrange[:] noextend
set yrange[:] noextend

# set size ratio -1
set view map
set datafile separator ','

set term gif animate delay 100 size 2000, 1000

set output "graphs/result_flat.gif"
set key center tmargin

do for [n = 0 : count : step] {
    filename = sprintf("data/main/%d.csv", n)
    plot filename with image title sprintf("Sorry %d", n)
}
set output
