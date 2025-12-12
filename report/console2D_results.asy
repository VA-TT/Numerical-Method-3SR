settings.outformat = "pdf";
size(6cm, keepAspect = true);
defaultpen(linewidth(1));

file console = input("console2D_graph.dat").line();
real[][] data = console;
data = transpose(data);
real[] ITER = data[0];
real[] DELTAX = data[1];
real[] DELTAY = data[2];
real[] DISPLACEMENT = data[3];
real[] FORCE = data[4];

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

// Macro: Draw fixed joint 
void groundCircle(pair A, real r)
{
    // Use a thin pen for these decorations
    pen p = linewidth(0.75);
    pen fillcolor = blue + p;
    pen bordercolor = black;
    pen tangentcolor = gray + p;
    real hatchSpacing=r;


    // tangent x postion on the right
    real xtan = A.x - r;

    // tangent y position
    real ymin = A.y - 3.0*r;
    real ymax = A.y + 3.0*r;

    //draw tangent
    draw((xtan, ymin)--(xtan, ymax), tangentcolor);

    // Hatch
    for(real y = ymin; y <= ymax; y += hatchSpacing) {
        draw( (xtan, y) -- (xtan - r * 0.5, y - r* 0.5), tangentcolor );
    }

        // Circle
    filldraw(circle(A, r), fillcolor, bordercolor);
}

//Axis 
Label i1 = Label("$\vec{i_1}$", position=EndPoint, align=S);
Label i2 = Label("$\vec{i_2}$", position=EndPoint, align=W);
draw((0,0)--(2,0),arrow=Arrow(size=6), L = i1);
draw((0,0)--(0,2),arrow=Arrow(size=6), L = i2);


// Coordinate of nodes
real a = 10;
pair A = (0,0);
pair B = (0, a);
pair C1 = (a, a);
real r = a / 40;

// Vẽ các thanh
Label b1 = Label("$\vec{b_1'} = \vec{x'}$", position=MidPoint, align=2SE);
Label b2 = Label("$\vec{b_2'}$", position=MidPoint, align=2S);
draw(A--C1, red+dashed);
draw(B--C1, red+dashed);




// Joints
groundCircle(A,r);
groundCircle(B,r);
filldraw(circle(A, r), blue, black);
filldraw(circle(B, r), blue, black);
filldraw(circle(C1, r), blue, black);

int lastIndex = ITER.length-1;
pair C2 = C1 + (DELTAX[lastIndex], DELTAY[lastIndex]);

draw(A--C2, red, L = b1, arrow = MidArrow(TeXHead));
draw(B--C2, red, L = b2, arrow = MidArrow(TeXHead));

//Force
pair F = (C2.x+2,C2.y-2);
Label f = Label("$\vec{F}$", position=EndPoint, align=E);
draw(C2--F, purple, arrow=Arrow(TeXHead), L = f);
draw(C2--(C2.x, F.y-1),dashed+gray);
pair p1 = (C2+F)/2;
pair p2 = (C2.x, C2.y-1);
draw(arc(C2, p1, p2, direction = CW), L = Label("$\theta$", position=MidPoint, align = S));

// Names (drawn with contrasting color so they are visible over the blue fill)
dot("$A$", A, 2*S, blue);
dot("$B$", B, 2*NE, blue);
dot("$C$", C1, 2*NE, blue);
dot("$C'$", C2, 2*NE, blue);