let m = 1.0; //mass
let m0 = 1.0; //mass
let G = 1.0e6; //grqvitqtionql constant

//kinematics
let x = 75.0; //mass
let y = 20.0; //stifness
let vx = -1.8 * y; //initial length opf the link
let vy = -1.5 * x; //initial length opf the link
let ax = 0.0;
let ay = 0.0;

let dt = 0.005;

let x0;
let y0;



function setup()
{
    createCanvas(500, 500);
    background(200, 25, 100);

    x0 = 0.5 * width;
    y0 = 0.5 * height;
}

function draw()
{
    // Euler time integration
    x += vx * dt;
    y += vy * dt;
    
    let d = sqrt(x * x + y * y);
    let ux = x / d;
    let uy = y / d;

    let a = (G * m0 * m / (d * d)) / m;
    ax = -a * ux;
    ay = -a * uy;

    vx += ax * dt;
    vy += ay * dt;

    stroke(0);
    fill(255, 255, 0);
    circle(x0, y0, 25)
    strokeWeight(2);
    line(x0, y0, x0 + x, y0 - y);
}

