import graph;
size(400,300);

// Domain length and analytic solution for the example
real L = 1;
real u_analytic(real x) { return 1.5*x - 0.25*x^2; }

// Read two-column numeric data from fem1D_data.txt
real[] xn, un;
string[] lines = readlines("fem1D_data.txt");
for (int i = 0; i < lines.length; ++i) {
  string s = trim(lines[i]);
  if (s == "") continue;
  string[] parts = split(s);
  if (parts.length >= 2) {
    xn.push(real(parts[0]));
    un.push(real(parts[1]));
  }
}

// Draw analytic curve
real xmin = 0, xmax = L;
path analytic = graph(u_analytic, xmin, xmax, 200);
draw(analytic, blue+1);

// Axes
xaxis(Label("$x$"), BottomTop, Ticks(Step=0.2));
yaxis(Label("$u(x)$"), LeftRight, Ticks(Step=0.25));

// Plot numeric points and connect them
for (int i = 0; i < xn.length; ++i) dot((xn[i], un[i]), red+4);
path numericPath = (xn[0], un[0]);
for (int i = 1; i < xn.length; ++i) numericPath = numericPath -- (xn[i], un[i]);
draw(numericPath, dashed+red);

// Legends
label("$\mathrm{Analytic:\ }u(x)=1.5x-0.25x^2$", (0.6, max(un)+0.15));
label("$\mathrm{Numeric\ nodes}$", (0.6, max(un)+0.02));

