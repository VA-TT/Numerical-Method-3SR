# fem1D_compare.gp - plot analytic and numeric FEM nodal values
set terminal pngcairo size 800,600 enhanced font 'Sans,12'
set output 'fem1D_compare.png'
set title 'fem1D: analytic vs numeric (EA=10, f=5, F=10)'
set xlabel 'x'
set ylabel 'u(x)'
set grid
set xrange [0:1]
set yrange [0:*]

# analytic function: u(x) = 1.5*x - 0.25*x^2
u(x) = 1.5*x - 0.25*x**2

# data file (two columns: x u)
datafile = 'fem1D_data.txt'

plot u(x) with lines lw 2 lc rgb 'blue' title 'analytic', \
     datafile using 1:2 with points pt 7 ps 1.5 lc rgb 'red' title 'numeric nodes'

unset output
