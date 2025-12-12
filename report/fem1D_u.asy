import graph;
settings.outformat = "pdf";
size(200,150,IgnoreAspect);

// Domain length and analytic solution
real L = 1;
real u_analytic(real x) { return 1.5*x - 0.25*x^2; }

// FEM 2 nodes
file twoNode = input("fem1D_2nodes.txt").line();
real[][] a = twoNode;
a = transpose(a);
real[] x1 = a[0];
real[] u1 = a[1];
path p1 = graph(x1,u1);
draw(p1, purple+dashed, "FEM $n=2$");


// FEM 6 nodes
file sixNode = input("fem1D_6nodes.txt").line();
real[][] b = sixNode;
b = transpose(b);
real[] x2 = b[0];
real[] u2 = b[1];
path p2 = graph(x2,u2);
draw(p2, green+dashed, "FEM $n=6$");
// Draw dots at FEM nodes
for(int i=0; i < x2.length; ++i){
    dot((x2[i], u2[i]), green);
}
for(int i=0; i < x1.length; ++i){
    dot((x1[i], u1[i]), purple);
}

// Analytic solution
path analytic = graph(u_analytic, 0, L, 200);
draw(analytic, red, "$1.5x - 0.25x^2$");

xaxis("$x$ (m)",BottomTop,LeftTicks);
yaxis("$u(x)$ (m)",LeftRight,RightTicks);

// Legend inside
add(scale(0.5)*legend(), (0.5,0.38));
