let softBody;
let floor;
let left;
let right;
let polygons = [];

let pause = false;
let mvx, mvy;
let mouseInteractionRadius = 0;

function setup() {
  createCanvas(800, 600);
  softBody = new Soft_body(350, 10, 15, 12, SPRING_CONSTRAINT_STIFFNESS, SOLVER_ITERATIONS);

  floor = new Polygon([createVector(0, height - 50), createVector(width, height - 50), createVector(width, height - 50 - 50), createVector(0, height - 50 - 50)]);
  left = new Polygon([createVector(0, 0), createVector(50, height), createVector(100, height), createVector(50, 0)]);
  right = new Polygon([createVector(width, 0), createVector(width - 50, 0), createVector(width - 50, height), createVector(width, height)]);
  ceiling = new Polygon([createVector(0, 0), createVector(width, 0), createVector(width, -50), createVector(0, -50)]);
  // polygons.push(new Polygon(500, 150, 9));
  // polygons.push(new Polygon(100, 300, 7));
  // polygons.push(new Polygon(350, 400, 5));
  polygons.push(ceiling);
  polygons.push(floor);
  polygons.push(left);
  polygons.push(right);
}

let dt = 10;

function draw() {
  // frameRate(0.2);
  background(51);
  push();
    noStroke();
    fill(201);
    textAlign(CENTER, TOP);
    textFont('Product Sans');
    textSize(15);    
    text("Click and drag with mouse to interact with the softbody", width / 2, 30);
  pop();


  if (mouseIsPressed && mouseButton === LEFT) {
    mvx = (mouseX - pmouseX) * MOUSE_FORCE_MULTPLIER;
    mvy = (mouseY - pmouseY) * MOUSE_FORCE_MULTPLIER;

  } else {
    mvx = 0;
    mvy = 0;
  }
  // dt = deltaTime
  
  // apply all external forces on all particles
  softBody.applyExternalForce(0, GRAVITY_FORCE_Y, dt);
  softBody.applyExternalForce(mvx, mvy, dt);
  softBody.dampVelocity(DAMPING);
  // noLoop();
  
  // calculate the predicted postions of all particles
  softBody.projectPosition(dt);
  
  for (let polygon of polygons){
    polygon.show();
  }
  
  softBody.generateCollisionConstraints(polygons);
  softBody.show();
  // solve step
  for (let i=0; i < SOLVER_ITERATIONS; i++) {
    // solve all constraints
    softBody.solveSpringConstraints();
    softBody.solveCollisionConstraints();
  }
  
  if (softBody.collisionConstraints.length > 0) {
    console.log(`deltaTime: ${deltaTime}, dt: ${dt}, constraints: ${softBody.collisionConstraints.length}`);
  }
  //post-solve step
  softBody.updateVelocity(dt);
  softBody.updatePosition();

  softBody.updateCollidingMassVelocity(RESTITUTION, FRICTION);
}
