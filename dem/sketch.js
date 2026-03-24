let x = 0;

function setup()
{
    createCanvas(500, 300);
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

