settings.outformat = "pdf";
size(10cm, keepAspect=true);
import graph;

// Read iteration vs deltaX
file itf = input("console2D_iterations.txt").line();
real[][] itA = itf; itA = transpose(itA);
real[] iter = itA[0];
real[] dx = itA[1];
real[] dy = itA[2];
real[] mag = itA[3];

// Create an axis picture for convergence
picture conv;
size(conv, 8cm, 4cm, keepAspect=true);
real maxIter = iter[iter.length-1];
real maxMag = max(mag);

draw(conv, xscale(1)*unitsquare, invisible); // initialize

// draw dx, dy, magnitude vs iteration
path p_dx = graph(iter, dx);
path p_dy = graph(iter, dy);
path p_mag = graph(iter, mag);
draw(conv, p_dx, blue+linewidth(1), "dx");
draw(conv, p_dy, green+linewidth(1), "dy");
draw(conv, p_mag, red+linewidth(1), "|d|\n");

// axes for the convergence plot
xaxis(conv, "Iterations", BottomTop, LeftTicks);
yaxis(conv, "Delta components / magnitude", LeftRight, RightTicks);

// Now read positions and draw structure before/after
file posf = input("console2D_positions.txt").line();
real[][] P = posf; P = transpose(P);
// P columns: name? Asymptote reads all tokens; assume order: name x0 y0 x1 y1 per line
// We'll read by value: x_before y_before x_after y_after

real Ax0 = P[0][0], Ay0 = P[1][0];
real Bx0 = P[0][1], By0 = P[1][1];
real Cx0 = P[0][2], Cy0 = P[1][2];
real Ax1 = P[2][0], Ay1 = P[3][0];
real Bx1 = P[2][1], By1 = P[3][1];
real Cx1 = P[2][2], Cy1 = P[3][2];

picture geom;
size(geom, 8cm, 6cm, keepAspect=true);
// draw before (dashed gray)
draw(geom, (Ax0,Ay0)--(Cx0,Cy0), dashed+gray);
draw(geom, (Bx0,By0)--(Cx0,Cy0), dashed+gray);
dot(geom, (Ax0,Ay0), blue);
dot(geom, (Bx0,By0), blue);
dot(geom, (Cx0,Cy0), blue);
label(geom, "before", (min(Ax0,Bx0,Cx0)+0.2, max(Ay0,By0,Cy0)+0.2));

// draw after (solid red)
draw(geom, (Ax1,Ay1)--(Cx1,Cy1), red+linewidth(1.5));
draw(geom, (Bx1,By1)--(Cx1,Cy1), red+linewidth(1.5));
dot(geom, (Ax1,Ay1), black);
dot(geom, (Bx1,By1), black);
dot(geom, (Cx1,Cy1), black);
label(geom, "after", (min(Ax1,Bx1,Cx1)+0.2, max(Ay1,By1,Cy1)+0.2));

// Compose final page
add(shift(0,3cm)*conv);
add(shift(0,-3cm)*geom);
