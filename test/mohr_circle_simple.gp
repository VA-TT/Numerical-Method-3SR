reset

# Usage:
#   gnuplot -e "sigma1=2; sigma3=0; c=1000; phi=35" test/mohr_circle_simple.gp
# Output:
#   mohr_circle_simple.pdf

sigma1 = 40000.0
sigma3 = 15000.0
c      = 1000.0
phi    = 35.0   # degrees


phi_rad = phi*pi/180.0
m = tan(phi_rad)

center = 0.5*(sigma1 + sigma3)
radius = 0.5*(sigma1 - sigma3)

set term pdfcairo enhanced color font 'Arial, 18' size 16cm,10cm
set output 'mohr_circle_simple.pdf'

set grid
set size ratio -1
unset key

set xlabel '{/Symbol s} (kPa)'
set ylabel '{/Symbol t} (kPa)'

# xmin = sigma3 - 0.2*radius
# xmax = sigma1 + 0.2*radius
# set xrange [xmin:xmax]

# ymin = -1.2*radius
# yline = c + m*xmax
# ymax = (1.2*radius > yline) ? (1.2*radius) : yline
# set yrange [ymin:ymax]

set object 1 circle at first center, 0 size radius fillstyle empty border lc rgb "#000000" lw 2

mc(x) = m*x + c
plot mc(x) with lines lc rgb "#e51e10" lw 2 notitle

set output
