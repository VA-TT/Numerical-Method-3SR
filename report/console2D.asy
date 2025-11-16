settings.outformat = "pdf";
size(6cm, keepAspect = true);
defaultpen(linewidth(1));

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

// Macro: Vẽ circle + tiếp tuyến dọc + nền gạch chéo
void groundCircle(pair A, real r)
{
    // Use a thin pen for these decorations
    pen p = linewidth(0.75);
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
    for(real y = ymin; y <= ymax; y += hatchSpacing) {
        draw( (xtan, y) -- (xtan - r * 0.5, y - r* 0.5), tangentcolor );
    }

        // Vẽ đường tròn
    filldraw(circle(A, r), fillcolor, bordercolor);
}

//Axis 
Label i1 = Label("$\vec{i_1}$", position=EndPoint, align=S);
Label i2 = Label("$\vec{i_2}$", position=EndPoint, align=W);
draw((0,0)--(2,0),arrow=Arrow(size=6), L = i1);
draw((0,0)--(0,2),arrow=Arrow(size=6), L = i2);


// Toạ độ các nút
real a = 10;
pair A = (0,0);
pair B = (0, a);
pair C = (a, a);
real r = a / 40;

// Vẽ các thanh
Label b1 = Label("$\vec{b_1} = \vec{x}$", position=MidPoint, align=2SE);
Label b2 = Label("$\vec{b_2}$", position=MidPoint, align=2S);
draw(A--C, red, L = b1, arrow = MidArrow(TeXHead));
draw(B--C, red, L = b2, arrow = MidArrow(TeXHead));

// Names (drawn with contrasting color so they are visible over the blue fill)
dot("$A$", A, 2*S, blue);
dot("$B$", B, 2*NE, blue);
dot("$C$", C, 2*NE, blue);

//Force
pair F = (C.x+2,C.y-2);
Label f = Label("$\vec{F}$", position=EndPoint, align=E);
draw(C--F, purple, arrow=Arrow(TeXHead), L = f);
draw(C--(C.x, F.y-1),dashed+gray);
pair p1 = (C+F)/2;
pair p2 = (C.x, C.y-1);

draw(arc(C, p1, p2, direction = CW), L = Label("$\theta$", position=MidPoint, align = S));

//Dim
pair dimShiftX = (-1cm, 0);
pair dimShiftY = (0, 1cm);
Label barLengthH = Label("$a$", align=N, position=MidPoint, filltype=Fill(white));
Label barLengthV = Label("$a$", align=W, position=MidPoint, filltype=Fill(white));
drawshifted(A -- B, trueshift=dimShiftX, label=barLengthV, arrow=Arrows(size=8), bar=Bars);
drawshifted(B -- C, trueshift=dimShiftY, label=barLengthH, arrow=Arrows(size=8), bar=Bars);

// Khớp (chấm nhỏ màu xanh)
groundCircle(A,r);
groundCircle(B,r);
filldraw(circle(A, r), blue, black);
filldraw(circle(B, r), blue, black);
filldraw(circle(C, r), blue, black);