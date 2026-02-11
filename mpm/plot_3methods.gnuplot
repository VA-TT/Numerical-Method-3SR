#!/usr/bin/gnuplot
# Comprehensive comparison of C++, Python, and Analytical solutions

set terminal pngcairo enhanced font 'Arial,11' size 1800,1200
set output 'comparison_3methods.png'

set multiplot layout 2,3 title "MPM 1D: C++ vs Python vs Analytical" font ",16"

# Plot 1: Position comparison
set title "Position vs Time"
set xlabel "Time (s)"
set ylabel "Position (m)"
set grid
set key right bottom
plot 'mpm1D_history.txt' using 1:4 with lines linewidth 3 title "Analytical", \
     'mpm1D_history.txt' using 1:2 with lines linewidth 2 dashtype 2 title "C++", \
     'mpm1D_python.txt' using 1:2 with points pointtype 7 pointsize 0.3 title "Python"

# Plot 2: Velocity comparison
set title "Velocity vs Time"
set xlabel "Time (s)"
set ylabel "Velocity (m/s)"
set grid
set key right top
plot 'mpm1D_history.txt' using 1:5 with lines linewidth 3 title "Analytical", \
     'mpm1D_history.txt' using 1:3 with lines linewidth 2 dashtype 2 title "C++", \
     'mpm1D_python.txt' using 1:3 with points pointtype 7 pointsize 0.3 title "Python"

# Plot 3: Phase space
set title "Phase Space (v vs x)"
set xlabel "Position (m)"
set ylabel "Velocity (m/s)"
set grid
set key right top
set size ratio -1
plot 'mpm1D_history.txt' using 4:5 with lines linewidth 3 title "Analytical", \
     'mpm1D_history.txt' using 2:3 with lines linewidth 2 dashtype 2 title "C++", \
     'mpm1D_python.txt' using 2:3 with points pointtype 7 pointsize 0.3 title "Python"
unset size

# Plot 4: Position error
set title "Position Error (vs Analytical)"
set xlabel "Time (s)"
set ylabel "Absolute Error (m)"
set grid
set logscale y
set key right top
plot 'mpm1D_history.txt' using 1:(abs($2-$4)) with lines linewidth 2 title "C++", \
     'mpm1D_python.txt' using 1:(abs($2-$4)) with lines linewidth 2 dashtype 2 title "Python"

# Plot 5: Velocity error
set title "Velocity Error (vs Analytical)"
set xlabel "Time (s)"
set ylabel "Absolute Error (m/s)"
set grid
set logscale y
set key right top
plot 'mpm1D_history.txt' using 1:(abs($3-$5)) with lines linewidth 2 title "C++", \
     'mpm1D_python.txt' using 1:(abs($3-$5)) with lines linewidth 2 dashtype 2 title "Python"

# Plot 6: Relative errors
set title "Relative Error (%)"
set xlabel "Time (s)"
set ylabel "Relative Error (%)"
set grid
unset logscale y
set key right top
set yrange [0:*]
plot 'mpm1D_history.txt' using 1:(abs($2-$4)/(abs($4)+1e-10)*100) with lines linewidth 2 title "C++ Position", \
     'mpm1D_python.txt' using 1:(abs($2-$4)/(abs($4)+1e-10)*100) with lines linewidth 2 dashtype 2 title "Python Position"

unset multiplot

print "3-method comparison plot generated: comparison_3methods.png"
