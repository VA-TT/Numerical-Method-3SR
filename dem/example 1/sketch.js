

let trajectory = [];
let t = 0;
let dt = 0.02;
let vx0 = 40.0;
let vy0 = 80.0;
let gx = 0;
let gy = -9.81;
let ballRadius = 8;
let simDone = false;

function run(vx0, vy0, dt) {
    let x = 0;
    let y = 0;
    let vx = vx0;
    let vy = vy0;
    trajectory = [];
    while (y >= 0.0 && x <= width) {
        trajectory.push({x: x, y: y});
        // Euler time integration
        x += vx * dt;
        y += vy * dt;
        vx += gx * dt;
        vy += gy * dt;
    }
    simDone = true;
}


function setup() {
    createCanvas(500, 300);
    background(240);
    run(vx0, vy0, dt);
    t = 0;
}

function draw() {
    background(255);

    // Vẽ quỹ đạo
    stroke(0, 0, 255);
    strokeWeight(2);
    noFill();
    beginShape();
    for (let i = 0; i < trajectory.length; i++) {
        vertex(trajectory[i].x, height - trajectory[i].y);
    }
    endShape();

    // Vẽ viên đạn di chuyển theo quỹ đạo
    if (trajectory.length > 0) {
        let idx = min(floor(t), trajectory.length - 1);
        let px = trajectory[idx].x;
        let py = trajectory[idx].y;
        fill(255, 0, 0);
        noStroke();
        ellipse(px, height - py, ballRadius * 2, ballRadius * 2);
        if (idx < trajectory.length - 1) {
            t += 1;
        }
    }
}

