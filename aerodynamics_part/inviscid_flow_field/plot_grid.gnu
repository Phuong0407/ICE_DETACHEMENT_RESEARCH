set title 'Computational Mesh'
set xlabel 'X-axis'
set ylabel 'Y-axis'
set grid
set key off
plot 'grid_points.dat' with points pt 7 ps 0.1 lc rgb 'blue', \
     'grid_connections.dat' with lines lc rgb 'black'
