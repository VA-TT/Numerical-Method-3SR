settings.outformat = "pdf";
size(8cm, keepAspect = true);

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

//Axis 
// Label i1 = Label("$\vec{i_1}$", position=EndPoint, align=S);
// Label i2 = Label("$\vec{i_2}$", position=EndPoint, align=W);
// draw((0,0)--(2,0),arrow=Arrow(size=6), L = i1);
// draw((0,0)--(0,2),arrow=Arrow(size=6), L = i2);


// Toạ độ các nút
real a = 10;
pair pA = (0,0);
pair pB = (a, 0);
real r = a / 50;

// Vẽ các thanh
draw(pA--pB, red, arrow = MidArrow(TeXHead));


// Concentrated force at node B (arrow pointing to the node)
pair F = (pB.x + 2, pB.y);
Label f = Label("$\vec{F}$", position=EndPoint, align=2SE);
// Draw concentrated force arrow at B (omit explicit margin to avoid name clash)
draw(F--pB, purple, arrow=Arrow(), L = f, margin=ArrowMargins);

//distributed borce
// draw(C--(C.x, F.y-1),dashed+gray);

// Dim
pair dimShiftX = (0.5cm, 0);
pair dimShiftY = (0, -1cm);

Label barLengthV = Label("$l$", align=N, position=MidPoint, filltype=Fill(white));
// use the correct parameter name 'label' (drawshifted's parameter is named 'label')
drawshifted(pA -- pB, trueshift=dimShiftY, label = barLengthV, arrow=Arrows(size=8), bar=Bars);

// Distributed load along the beam from A to B (spacing = 2)
distributedF(pA, pB, 1);
dot(pB, blue);

// Khớp (chấm nhỏ màu xanh) - supports at A and B
groundCircleV(pA, r);

// Labels for nodes
label("$A$", pA, 2*NE);
label("$B$", pB, 2*NE);
