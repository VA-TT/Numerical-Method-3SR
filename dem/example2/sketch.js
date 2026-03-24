let m = 0.005; //mass
let k = 150.0; //stifness
let d0 = 100.0; //initial length opf the link
let gy = 4000.0;
let v0 = 100.0;


let x = 0.0; //mass
let y = d0; //stifness
let vx = v0; //initial length opf the link
let ax = 0.0;
let ay = 0.0;

let dt = 0.005;

let x0;
let y0;

function run(vx0, vy0, dt)
{
    let x = 0;
    let y = 0;
    let vx = vx0;
    let vy = vy0;
    let gx = 0;
    let gy = -9.81;

    stroke(0);
    strokeWeight(1.5);

    stroke(255, 0, 0);
    strokeWeight(2);
    line(0, height, vx0, height - vy0);
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


        let prev_x = x;
        let prev_y = y;

        //Euler time integration
    x += vx * dt;
    y += vy * dt;
    
    let d = sqrt(x * x + y * y);
    ux = x / d;
    uy = y / d;

    let a = (k * (d - d0) / m);
    ax = -a * ux;
    ay = -a * uy - gy;

    vx += gx * dt;
    vy += gy * dt;

    //damp
    vx *= 0.995;
    vy *= 0.995;
    
    fill(236, 208, 247);
    stroke(0);
    strokeWeight(2);
    line(x0, y0, x0 + x, y0 - y);
}

