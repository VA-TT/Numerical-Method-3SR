
let m = 0.005; // mass
let k = 150.0; // stiffness
let d0 = 100.0; // initial length of the link
let gy = 4000.0; // gravity acceleration
let v0 = 100.0; // initial velocity

let x = 0.0; // position x
let y = d0; // position y
let vx = v0; // velocity x
let vy = 0.0; // velocity y (added)
let ax = 0.0;
let ay = 0.0;

let dt = 0.005;

let x0;
let y0;

function run(vx0, vy0, dt)
{
    draw();
}


function setup()
{
    createCanvas(500, 500);
    background(200, 25, 100);

    x0 = 0.5 * width;
    y0 = 0.5 * height;
}

function draw()
{
    background(255);

    // Euler time integration
    x += vx * dt;
    y += vy * dt;

    let d = sqrt(x * x + y * y);
    let ux = x / d;
    let uy = y / d;

    let a = (k * (d - d0) / m);
    ax = -a * ux;
    ay = -a * uy - gy;

    vx += ax * dt;
    vy += ay * dt;

    // Damping
    vx *= 0.995;
    vy *= 0.995;

    fill(236, 208, 247);
    stroke(0);
    strokeWeight(2);
    line(x0, y0, x0 + x, y0 - y);
}

