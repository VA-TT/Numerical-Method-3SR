set terminal pngcairo size 1200,900
set output './gradient_scalar_z0.png'

set title 'Contour + Gradient field of f(x,y,z)=1.5xy+2z on z=0'
set xlabel 'x'
set ylabel 'y'
set xrange [0:1]
set yrange [0:1]
set size ratio -1
set grid
set key off

# Build contour lines from scalar field f in column 6.
set view map
set contour base
unset surface
set cntrparam levels 12
set table $contours
splot './gradient_scalar_z0.dat' using 1:2:6
unset table

# Restore 2D plotting and style settings.
set palette rgb 33,13,10
set cblabel '|grad_{xy}|'

# Scale arrow lengths for readability on [0,1]x[0,1]
s = 0.18

# Data columns: x y gx gy |grad_xy| f
plot \
	$contours using 1:2 with lines lw 1.2 lc rgb '#202020', \
	'./gradient_scalar_z0.dat' using 1:2:($3*s):($4*s):5 with vectors filled head lw 1.2 lc palette
