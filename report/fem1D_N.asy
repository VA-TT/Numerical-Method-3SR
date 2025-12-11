import graph;
settings.outformat = "pdf";
size(200,150,IgnoreAspect);

// Domain length and analytic solution
real L = 1;
real N_analytic(real x) { return 15 - 5*x; }

// FEM 2 nodes
file twoNode = input("fem1D_2nodes.txt").line();
real[][] a = twoNode;
a = transpose(a);
real[] x1 = a[0];
real[] u1 = a[1];
real[] N1 = a[2];
path p1 = graph(x1,N1);
draw(p1, purple+dashed, "FEM $n=2$");
for(int i=0; i < x1.length; ++i){
    dot((x1[i], N1[i]), purple);
}

// FEM 5 nodes
file sixNode = input("fem1D_6nodes.txt").line();
real[][] b = sixNode;
b = transpose(b);
real[] x2 = b[0];
real[] u2 = b[1];
real[] N2 = b[2];
path p2 = graph(x2,N2);
draw(p2, green+dashed, "FEM $n=6$");
// Draw dots at FEM nodes
for(int i=0; i < x2.length; ++i){
    dot((x2[i], N2[i]), green);
}

// Analytic solution
path analytic = graph(N_analytic, 0, L, 200);
draw(analytic, red, "$15 - 5x$");

xaxis("$x$ (m)",BottomTop,LeftTicks);
yaxis("$N(x)$ (kN)",LeftRight,RightTicks);

add(scale(0.5)*legend(), (0.1,11.5));