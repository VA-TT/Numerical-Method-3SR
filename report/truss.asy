settings.outformat = "pdf";
size(8cm, keepAspect = true);

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

// Macro: Vẽ circle + tiếp tuyến ngang + nền gạch chéo
void groundCircle(pair A, real r)
{
    // Use a thin pen for these decorations
    pen p = linewidth(0.75);
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
    for(real x = xmin; x <= xmax; x += hatchSpacing) {
        draw( (x, ytan) -- (x - r* 0.5, ytan - r * 0.5), tangentcolor );
    }

        // Vẽ đường tròn
    filldraw(circle(A, r), fillcolor, bordercolor);
}

//Axis 
// Label i1 = Label("$\vec{i_1}$", position=EndPoint, align=S);
// Label i2 = Label("$\vec{i_2}$", position=EndPoint, align=W);
// draw((0,0)--(2,0),arrow=Arrow(size=6), L = i1);
// draw((0,0)--(0,2),arrow=Arrow(size=6), L = i2);


// Toạ độ các nút
real a = 10;
pair pA = (-a,0);
pair pB = (a, 0);
pair pC = (-a, 2a);
pair pD = (a, 2a);
pair pE = (-a, 4a);
pair pF = (a, 4a);
real r = a / 10;

// Vẽ các thanh
draw(pA--pC, red, arrow = MidArrow(TeXHead), L = Label("$\vec{b_1}$", position=MidPoint, align=W));
draw(pC--pB, red, arrow = MidArrow(TeXHead), L = Label("$\vec{b_2}$", position=MidPoint, align=2S));
draw(pB--pD, red, arrow = MidArrow(TeXHead), L = Label("$\vec{b_3}$", position=MidPoint, align=W));
draw(pC--pD, red, arrow = MidArrow(TeXHead), L = Label("$\vec{b_4}$", position=MidPoint, align=N));
draw(pC--pE, red, arrow = MidArrow(TeXHead), L = Label("$\vec{b_5}$", position=MidPoint, align=W));
draw(pD--pE, red, arrow = MidArrow(TeXHead), L = Label("$\vec{b_6}$", position=MidPoint, align=2S));
draw(pD--pF, red, arrow = MidArrow(TeXHead), L = Label("$\vec{b_7}$", position=MidPoint, align=W));
draw(pE--pF, red, arrow = MidArrow(TeXHead), L = Label("$\vec{b_8}$", position=MidPoint, align=N));


// //Force
// pair F = (C.x+2,C.y-2);
// Label f = Label("$\vec{F}$", position=EndPoint, align=E);
// draw(C--F, purple, arrow=Arrow(TeXHead), L = f);
// draw(C--(C.x, F.y-1),dashed+gray);
// pair p1 = (C+F)/2;
// pair p2 = (C.x, C.y-1);

// draw(arc(C, p1, p2, direction = CW), L = Label("$\theta$", position=MidPoint, align = S));

//Dim
// pair dimShiftX = (-1cm, 0);
// pair dimShiftY = (0, 1cm);
// Label barLengthH = Label("$a$", align=N, position=MidPoint, filltype=Fill(white));
// Label barLengthV = Label("$a$", align=W, position=MidPoint, filltype=Fill(white));
// drawshifted(A -- B, trueshift=dimShiftX, label=barLengthV, arrow=Arrows(size=8), bar=Bars);
// drawshifted(B -- C, trueshift=dimShiftY, label=barLengthH, arrow=Arrows(size=8), bar=Bars);

// Khớp (chấm nhỏ màu xanh)
groundCircle(pA,r);
groundCircle(pB,r);
filldraw(circle(pC, r), blue, black);
filldraw(circle(pD, r), blue, black);
filldraw(circle(pE, r), blue, black);
filldraw(circle(pF, r), blue, black);

// Sau đó vẽ label thì chúng sẽ nằm phía trên
label("$A$", pA, 2*NW);
label("$B$", pB, 2*NE);
label("$C$", pC, 2*NW);
label("$D$", pD, 2*NE);
label("$E$", pE, 2*NW);
label("$F$", pF, 2*NE);;