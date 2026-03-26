class Grain 
{
    constructor(mass, radius, xpos, ypos)
    {
        this.mass = mass;
        this.inertia = 0.5 * mass * radius * radius;
        this.radius = radius;

        this.xpos = xpos;
this.ypos = ypos;
        this.rot = 0.0;
        this.vx = 0.0;
        this.vy = 0.0;
        this.vrot = 0.0;
        this.fx = 0.0;
        this.fy = 0.0;
        this.frot = 0.0;

    }

    draw()
    {
        let vel = sqrt(this.vx * this.vx + this.vy * this.vy);
        let M = max(vRand, 0.8 * abs(topv));
        let colorVel = map(vel, 0.0, M, 255, 0);
        fill(colorVel);

        push();
		translate(this.xpos * scale, height - this.ypos * scale);
		rotate(-this.rot);
		stroke(0);
		let r = this.radius * scale;
		circle (0,0, 2*r);
		// line(0,0,r,0);
		pop();
        
    }
    }
	
let grain = [];
let ngw = 5;
let ngh = 5;
let rmin = 0.5e-3;
let rmax = 1e-3;
let density = 2700.0;
let scale = 40000;
let xmax = 0.0;
let ymax = 0.0;
let ymax0 = 0.0;
let vRand; 

let ng = ngh * ngw;
let meanDiameter = (rmin + rmax)/2.0;
let meanMass = density * 3.1415 * meanDiameter * meanDiameter / 4.0;

//kinemactic of the wall
let topf = 0.0;
let topa = 0.0;
let topv = 0.0;
let rightf = 0.0;
let righta = 0.0;
let rightv = 0.0;

//neighbor array
let neighbor = [];
let dmax = 0.95 * rmin;
let nStepsUpdate = 50;

//Processing
let fnmax = 0.0;
let viscoRate = 0.95;

let kn = 10000.0;
let kt = kn;
let mu = 0.5;
let fcohesion = 0.0;

let iStep = 0;

let gx = 0.0;
let gy = -9.81;

let dt;

function putOnGrid()
{
	let Vsolid = 0.0;
	for (let iy = 0; iy < ngh; ++iy)
	{
			for (let ix = 0; ix < ngw; ++ix)
	{
		let x = rmax + ix * 2 * rmax;
		let y = rmax + iy * 2 * rmax;
		let radius = random(rmin, rmax);
		let surf = PI * radius * radius;
		Vsolid += surf;
		let mass = surf * density;
		grain.push(new Grain(mass, radius, x ,y));
	}
	}
	dt = 3.14159 * sqrt(meanMass / kn) / 25.0;
	vRand = rmin / (500.0 * dt);
	for (let i = 0; i < grain.length; ++i)
	{
		grain[i].vx = random(-vRand, vRand);
		grain[i].vy = random(-vRand, vRand);
	}
}



function setup()
{
    createCanvas(500, 500);
    background(240);
	putOnGrid();
	// updateNeighborList();
	// noLoop();
}

function contact(){}

function computeForceBottomWall(i)
{ 

	let dn = grain[i].ypos - grain[i].radius;
	if(dn < 0.0)	{	
		let vn = -grain[i].vx;
		let visco = viscoRate * 2.0 * sqrt(grain[i].mass * kn);
		let fn = -kn * dn - visco*vn;
		grain[i].fy += fn;
		// grain[i].p += fn;
		
	}
}
function periodicBoundary(){
}

function updateAcceleration()
{
		for (let i = 0; i < grain.length; ++i)
	{
		grain[i].ax = gx;
		grain[i].ay = gy;
}
}
function timeIntegration(){
	updateAcceleration();

	for (let i = 0; i < grain.length; ++i)
	{	
		grain[i].xpos += grain[i].vx * dt;
		grain[i].ypos += grain[i].vy * dt;
		grain[i].vx += grain[i].ax * dt;
		grain[i].vy += grain[i].ay * dt;
		computeForceBottomWall(i);
	}
	
		for (let i = 0; i < grain.length; ++i)
	{	
		grain[i].xpos += grain[i].vx * dt;
		grain[i].ypos += grain[i].vy * dt;
		grain[i].vx += grain[i].ax * dt;
		grain[i].vy += grain[i].ay * dt;
		computeForceBottomWall(i);
	}

}

function draw()
{	background(240)
    stroke(255, 0, 0);
    strokeWeight(2);
	for (let i = 0; i < grain.length; ++i)
	{
		grain[i].draw();
	}
	timeIntegration();
}