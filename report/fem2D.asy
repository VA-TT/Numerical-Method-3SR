settings.outformat = "pdf";
size(8cm, keepAspect = true);

defaultpen(fontsize(8pt));

real lw = linewidth(currentpen);
// Margin object to shorten paths so arrowheads don't overlap endpoints
margin ArrowMargins = TrueMargin(lw, 0.2*lw);

//Macro draw shifted
void drawshifted(path g, pair trueshift, picture pic =
currentpicture, Label label="", pen pen=currentpen,
arrowbar arrow=None, arrowbar bar=None, margin
margin=NoMargin, marker marker=nomarker)
{
pic.add(new void(frame f, transform t) {
picture opic;
draw(opic, L=label, shift(trueshift)*t*g, p=pen,
arrow=arrow, bar=bar,
margin=margin, marker=marker);
add(f,opic.fit());
});
pic.addBox(min(g), max(g), trueshift+min(pen),
trueshift+max(pen));
}

// Macro: Vẽ circle + tiếp tuyến ngang + nền gạch chéo
void groundCircleH(pair A, real r)
{
    // Use a thin pen for these decorations
    pen p = linewidth(3/4 * lw);
    pen fillcolor = blue + p;
    pen bordercolor = black;
    pen tangentcolor = gray + p;
    real hatchSpacing=r;


    // Tính tọa độ tiếp tuyến dọc bên phải
    real ytan = A.y - r;

    // Chiều rong hiển thị cho tiếp tuyến và hatch
    real xmin = A.x - 3.0*r;
    real xmax = A.x + 3.0*r;

    // Vẽ đường tiếp tuyến màu ghi
    draw((xmin, ytan)--(xmax, ytan), tangentcolor);

    // Vẽ các đường hatch chéo giống ký hiệu đất
    for(real x = xmin+0.5*hatchSpacing; x <= xmax-0.5*hatchSpacing; x += hatchSpacing) {
        draw( (x, ytan) -- (x - r* 0.5, ytan - r * 0.5), tangentcolor );
    }

        // Vẽ đường tròn
    filldraw(circle(A, r), fillcolor, bordercolor);
}


// Macro: Vẽ circle + tiếp tuyến dọc + nền gạch chéo
void groundCircleV(pair A, real r)
{
    // Use a thin pen for these decorations
    pen p = linewidth(3/4 * lw);
    pen fillcolor = blue + p;
    pen bordercolor = black;
    pen tangentcolor = gray + p;
    real hatchSpacing=r;


    // Tính tọa độ tiếp tuyến dọc bên phải
    real xtan = A.x - r;

    // Chiều cao hiển thị cho tiếp tuyến và hatch
    real ymin = A.y - 3.0*r;
    real ymax = A.y + 3.0*r;

    // Vẽ đường tiếp tuyến màu ghi
    draw((xtan, ymin)--(xtan, ymax), tangentcolor);

    // Vẽ các đường hatch chéo giống ký hiệu đất
    for(real y = ymin+0.5*hatchSpacing; y <= ymax-0.5*hatchSpacing; y += hatchSpacing) {
        draw( (xtan, y) -- (xtan - r * 0.5, y - r* 0.5), tangentcolor );
    }

        // Vẽ đường tròn
    filldraw(circle(A, r), fillcolor, bordercolor);
}

//Macro draw distributeed force
void distributedF(pair pointStart, pair pointEnd, real f)
{
    // Use a thin pen for these decorations
    pen p = linewidth(3/4 * lw);
    pen forceColor = purple;
    real hatchSpacing= abs((pointStart.x - pointEnd.x))/5;


    // Boundary -- use Label constructor (Label(...)) for the draw L parameter
    draw((pointStart.x, pointStart.y+f)--(pointEnd.x, pointEnd.y+f), forceColor, L = Label("$\vec{f}$", position = MidPoint, align = N));

    // arrow forces
    for(real x = pointStart.x; x <= pointEnd.x; x += hatchSpacing) {
           draw( (x, pointStart.y +f) -- (x, pointStart.y), forceColor, arrow = Arrow(TeXHead), margin = ArrowMargins );
    }
}




// Toạ độ các nút
real a = 10;
pair pA = (10,0);
pair pB = (5, 10);
pair pC = (0, 10);
pair pO = (0, 0);

label("$\Omega$", (pA+pB+pC+pO)/4);

// Vẽ các side
path object = pO--pA--pB--pC--cycle;
draw(object, red);



// Hàm vẽ pháp tuyến tại điểm P trên đoạn MN
void drawNormal(pair M, pair N) {
    pair P = (M + N)/2;              // điểm giữa đoạn MN
    pair v = N - M;                  // vector cạnh
    pair n = rotate(90)*v;           // vector pháp tuyến (vuông góc)
    n = n/abs(n);                    // chuẩn hóa

    draw(P -- (P - n), purple, arrow = Arrow(TeXHead), L = Label("$\vec{n}$", position = EndPoint));       // vẽ pháp tuyến dài 2 đơn vị
}
// // Draw normal vector
drawNormal(pA, pB);
drawNormal(pB, pC);
drawNormal(pC, pO);
drawNormal(pO, pA);

// ground hatch
real xmin = -10;
real xmax= 15;
real hatchSpacing = 1;
    for(real x = xmin+0.5*hatchSpacing; x <= xmax-0.5*hatchSpacing; x += hatchSpacing) {
        draw( (x, 0) -- (x - 1* 0.5, 0 - 1 * 0.5), gray);
    }
draw((xmin,0) -- (xmax,0), gray);
//Draw water
draw((xmin, pC.y-1) -- (pC.x, pC.y-1), cyan);

// Labels for nodes
label("$O$", pO, 2*NE);
label("$A$", pA, 2*NE);
label("$B$", pB, 2*NE);
label("$C$", pC, 2*NE);

// Axis 
Label i1 = Label("$\vec{i_1}$", position=EndPoint, align=S);
Label i2 = Label("$\vec{i_2}$", position=EndPoint, align=W);
draw((0,0)--(2,0),arrow=Arrow(size=6), L = i1);
draw((0,0)--(0,2),arrow=Arrow(size=6), L = i2);

label("$\Gamma_u$", (pO+pA)*2/3, N);
label("$\Gamma_F = AB \cup BC \cup CO$", midpoint(pO+pC)/4*3, 3*W);
