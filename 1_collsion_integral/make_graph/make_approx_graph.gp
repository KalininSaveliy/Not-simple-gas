# set terminal png
name = "data/relax_time.txt"
set output "graph/approx_time.png"
set title "Relaxation time"

set xlabel "t, время свободного пробега"
shift =  -2.28009
set ylabel sprintf("log(M_{2} %.2f)", shift)

# set xrange[0:20]
# set logscale y
# set yrange [1:5]
# set format y '%.1e'

set xtics 5
set mxtics 5
set ytics 0.5
# set mytics 1
set grid xtics ytics mxtics mytics


# fitting with polynom
set label 1 "y = k * x + b" at 6, -10
f1(x) = k1 * x + b1
fit [0:15] f1(x) name u 1:(log($2 + shift)) via k1, b1
set label 2 sprintf("k_{xx} = %f", k1) at 2, -11
set label 3 sprintf("b_{xx} = %f", b1) at 2, -12

f2(x) = k2 * x + b2
fit [0:15] f2(x) name u 1:(log($3 + shift)) via k2, b2

set label 5 sprintf("k_{yy} = %f", k2) at 9, -11
set label 6 sprintf("b_{yy} = %f", b2) at 9, -12
# end of fitting

plot name u 1:(log($2 + shift)) title sprintf("log(M_{xx} %.2f)", shift) w l lc rgb "red",\
     name u 1:(log($3 + shift)) title sprintf("log(M_{yy} %.2f)", shift) w l lc rgb "green",\
     name u 1:(log($4 + shift)) title sprintf("log(M_{zz} %.2f)", shift) w l lc rgb "blue",\
     f1(x) title "approx XX",\
     f2(x) title "approx YY and ZZ"
