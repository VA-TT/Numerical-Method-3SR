reset
set autoscale
unset log
unset label
set border 0

# ============================
# Cercle de Mohr (simple)
# Inputs are in kPa (edit here, or override with: gnuplot -e "sigma1=...; sigma3=...; c=...; phi=..." ...)
# ============================

sigma1 = 7.06   # kPa
sigma3 = 2.872   # kPa
c      = 1.0       # kPa
phi    = 35.0      # degrees

# Plot ranges (kPa)
xmin = 0
xmax = sigma1*1.10
ymin = 0
ymax = 0.5*(sigma1 - sigma3)*1.20

set term pdfcairo enhanced color font 'Arial, 18' size 16cm,10cm
set output 'mohr_circle_simple.pdf'

set arrow 1 from 0,0 to xmax,0 nohead lw 2
set arrow 2 from 0,0 to 0,ymax nohead lw 2

set xlabel "{/Symbol s}_n (kPa)"
set ylabel "{/Symbol t} (kPa)"
unset key
unset grid

set xtics nomirror
set ytics nomirror

set xrange [xmin:xmax]
set yrange [ymin:ymax]
set size ratio -1

center = (sigma1 + sigma3)/2.0
radius = (sigma1 - sigma3)/2.0

phi_rad = phi*pi/180.0
m = tan(phi_rad)
mc(x) = m*x + c

set object 1 circle at first center, 0 size radius fillstyle empty border lc rgb "#000000" lw 2
plot mc(x) with lines lc rgb "#e51e10" lw 2 notitle

set output
