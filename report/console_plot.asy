import graph;
settings.outformat = "pdf";
size(10cm,6cm,IgnoreAspect);
defaultpen(linewidth(1));

real msize = 2.5;

// Marker fill
marker[] MarkFill = {
  marker(scale(msize)*unitcircle,   blue, Fill),
  marker(scale(msize)*polygon(3),   Fill),
  marker(scale(msize)*polygon(4),   red, Fill),
  marker(scale(msize)*polygon(5),   Fill),
  marker(scale(msize)*(invert*polygon(3)), Fill),
  marker(scale(msize)*diamond,      Fill)
};

// Đọc dữ liệu
file console = input("console2D_graph.dat").line();
real[][] data = console;
data = transpose(data);
real[] ITER = data[0];
real[] DISPLACEMENT = data[3];
real[] FORCE = data[4];

// ================================
//   LOG SCALE CHO TRỤC CHÍNH
// ================================
scale(Linear, Log);   // X = tuyến tính, Y1 = log

// Vẽ displacement (Y1)
path p_disp = graph(ITER, DISPLACEMENT);
draw(p_disp, blue, "Displacement", marker=MarkFill[0]);

yaxis("Displacement (m)", Left, blue,
      RightTicks(n=4));
xaxis("Iteration", BottomTop, LeftTicks(n=4));

// ================================
//   LOG SCALE CHO TRỤC Y THỨ HAI
// ================================
picture q = secondaryY(new void(picture pic) {

    scale(pic, Linear, Log);   // X linear, Y2 log

    path p_force = graph(pic, ITER, FORCE);
    draw(pic, p_force, red, "Force", marker=MarkFill[2]);

    yaxis(pic, "Force (Pa)", Right, red,
          LeftTicks(n=3));
});

add(q);

// // ================================
// //   Legend
// // ================================
// add(scale(0.6)*legend(), (2.2,1e-7));
