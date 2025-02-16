# example of call
# gnuplot -c .\make_graph\one_time_flat.gp 20
filename = "data/main/".ARG1.".csv"
print "filename     : ", filename

# set terminal png

set output "graphs/result_one_time_flat.png"
set datafile separator ','

set xrange[:] noextend
set yrange[:] noextend
# set autoscale xfix
# set autoscale yfix
# set autoscale cbfix
set grid

set title "iter = ".ARG1
plot filename with image notitle

pause -0.1
