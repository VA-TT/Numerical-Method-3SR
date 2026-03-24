

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

    while (y >= 0.0 && x <= width)
    {
        let prev_x = x;
        let prev_y = y;

        //Euler time integration
        x += vx * dt;
        y += vy * dt;

        vx += gx * dt;
        vy += gy * dt;

        //draw line of trajectory
        line(prev_x, height - prev_y, x, height - y);
    }

    stroke(255, 0, 0);
    strokeWeight(2);
    line(0, height, vx0, height - vy0);
}


function setup()
{
    createCanvas(500, 300);
    background(240);

    run(20.0, 80.0, 40.0);
}

function draw()
{
    background(255);

    x = x + 10;
    if (x > 250) x = -250;

    fill(255, 0, 0);
    noStroke();
    circle(250 + x, 150, 110);
}

