import three;
import math;

settings.outformat = "pdf";
size(10cm, keepAspect = true);
settings.render = 8;

// ======================
// Camera
// ======================
currentprojection = perspective(
    camera = (22, -33, 17),
    up = Z
);
pen supportPen = gray(0.4) + linewidth(0.8);

void fixedSupport(triple P, real s=0.8) {
    triple dx = (s,0,0);
    triple dy = (0,s,0);

    draw(P - dx - dy -- P + dx + dy, supportPen);
    draw(P - dx + dy -- P + dx - dy, supportPen);
}

// ======================
// Node coordinates
// ======================
triple[] nodes = {
    (0, 0, 0),        // 0
    (10, 0, 0),       // 1
    (0, 0, 10),       // 2
    (10, 0, 10),      // 3
    (0, -10, 0),      // 4
    (10, -10, 0),     // 5
    (0, -10, 10),     // 6
    (10, -10, 10)     // 7
};

// Các nút ngàm
int[] fixedNodes = {0,1,4,5};

for (int i=0; i<fixedNodes.length; ++i) {
    fixedSupport(nodes[fixedNodes[i]]);
}

triple[] nodesDeformed = {
    (0, 0, 0),        
    (10, 0, 0),       
    (0, 0, 9.38083),  
    (10, 0, 9.38083), 
    (0, -10, 0),      
    (10, -10, 0),     
    (0, -10, 9.38083),
    (10, -10, 9.38083)
};

// ======================
// Bar connectivity
// ======================
int[][] bars = {
    {0,2}, {2,3}, {3,1},
    {2,6}, {3,7}, {4,6},
    {6,7}, {7,5}, {3,6}, {2,7}
};

// ======================
// Pens & arrows
// ======================
pen barPen   = blue + linewidth(1.2);
pen forcePen = purple + linewidth(1.5);
pen nodePen  = red;

arrowbar3 barArrow = Arrow3(TeXHead2, position=0.5);
     // vector thanh
arrowbar3 forceArrow = Arrow3(TeXHead2);      // lực

// ======================
// Draw bars as vectors
// ======================
for (int i = 0; i < bars.length; ++i) {
    int i0 = bars[i][0];
    int i1 = bars[i][1];

    triple A = nodes[i0];
    triple B = nodes[i1];
    triple C = nodesDeformed[i0];
    triple D = nodesDeformed[i1];

    // Vẽ vector thanh
    draw(C -- D, barPen, barArrow);
    draw(A -- B, barPen+dashed);

    // Nhãn vector
    triple mid = 0.5*(A + B);
    triple offset = (0, 1, 0);
    label("$\vec{b_" + string(i) + "}$", mid + offset);
}

// ======================
// Draw forces
// ======================
triple F = (0, 0, 3);   // hướng xuống
int[] loadNodes = {2, 3, 6, 7};

for (int i = 0; i < loadNodes.length; ++i) {
    int n = loadNodes[i];
    draw((nodes[n] + F -- nodes[n]),
         forcePen, forceArrow);
    label("$F$", nodes[n] + F, N);
}

// ======================
// Draw nodes
// ======================
for (int i = 0; i < nodesDeformed.length; ++i) {
    dot(nodes[i], nodePen + 6bp);
    dot(nodesDeformed[i], nodePen + 6bp);
    label(string(i), nodesDeformed[i], NE);
}

// ======================
// Axes
// ======================
draw((0,0,0)--(12,0,0), Arrow3);
draw((0,0,0)--(0,-12,0), Arrow3);
draw((0,0,0)--(0,0,15), Arrow3);

label("$x$", (12,0,0), E);
label("$y$", (0,-12,0), S);
label("$z$", (0,0,15), W);
