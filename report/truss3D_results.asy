settings.outformat = "png";
settings.render = 16;

import three;      // Phải có dòng này
import graph3;

size(12cm);
currentprojection = orthographic(4, -2, 3);

// Đọc file dữ liệu nút
file trussFile = input("truss3D_node.dat").line();
real[][] nodeData = trussFile;
nodeData = transpose(nodeData);

int nNodes = nodeData[0].length;
triple[] nodes;
for(int i = 0; i < nNodes; ++i) {
    nodes.push((nodeData[1][i], nodeData[2][i], nodeData[3][i]));
}

// Định nghĩa các thanh
int[] barOrigin = {0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7};
int[] barEnd    = {1, 2, 3, 2, 3, 0, 1, 0, 5, 6, 7, 6, 7, 4, 5, 4};
int nBars = barOrigin.length;

// Macro vẽ gối tựa 3D
void groundSphere(triple P, real r) {
    draw(circle(P, r), blue);
    
    // Mặt phẳng đất
    real dx = 2*r;
    triple A = P + (-dx, -dx, -r);
    triple B = P + ( dx, -dx, -r);
    triple C = P + ( dx,  dx, -r);
    triple D = P + (-dx,  dx, -r);
    
    draw(surface(A--B--C--D--cycle), gray(0.7)+opacity(0.5));
}

// Vẽ trục tọa độ
real axisLen = 12;
draw(Label("$x$", 1), O -- axisLen*X, red, Arrow3);
draw(Label("$y$", 1), O -- axisLen*Y, green+dashed, Arrow3);
draw(Label("$z$", 1), O -- axisLen*Z, blue, Arrow3);

// Vẽ các thanh
for(int b = 0; b < nBars; ++b) {
    draw(nodes[barOrigin[b]] -- nodes[barEnd[b]], red+linewidth(1.5pt));
}

// Vẽ nút
real nodeR = 0.4;
int[] fixedNodes = {0, 1, 4, 5};
bool[] isFixed = array(nNodes, false);
for(int i : fixedNodes) isFixed[i] = true;

for(int n = 0; n < nNodes; ++n) {
    if(isFixed[n]) {
        groundSphere(nodes[n], nodeR);
    } else {
        draw(circle(nodes[n], nodeR), cyan);
    }
    label("$" + string(n) + "$", nodes[n], 2*Z);
}

// Vẽ lực (F = -30kN theo z tại nút 2,3,6,7)
triple forceVec = (0, 0, -2.5);
int[] loadNodes = {2, 3, 6, 7};
for(int n : loadNodes) {
    draw(nodes[n] -- nodes[n] + forceVec, purple+linewidth(2pt), Arrow3);
    label("$\vec{F}$", nodes[n] + forceVec, S);
}