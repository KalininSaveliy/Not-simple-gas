count = 100
step = 10
max_z = 0.01

set xrange [-5:10]
set yrange [-5:10]

# set palette defined (0 "white", max_z "red")
set cbrange [0:max_z]

set size ratio -1
set view map

set term gif animate delay 100 size 2000, 1000

set output "graphs/result_flat.gif"

do for [n = 0 : count : step] {
    # filename = sprintf("data/task3/out_%03d.dat", n)
    filename = sprintf("data/main/%d.csvs", n)
    filename = "data/main/0.csv"
    splot filename with image title sprintf("Sorry %d", n)
}
