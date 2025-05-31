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
  softBody = new Soft_body(350, 10, 10, 7, SPRINT_CONSTRAINT_STIFFNESS, SOLVER_ITERATIONS);
  floor = new Polygon([createVector(0, height - 50), createVector(width, height - 50), createVector(width, height - 50 - 50), createVector(0, height - 50 - 50)]);
  left = new Polygon([createVector(0, 0), createVector(50, height), createVector(100, height), createVector(50, 0)]);
  right = new Polygon([createVector(width - 150, 0), createVector(width - 150, height), createVector(width-80, height), createVector(width-80, 0)]);
  polygons.push(new Polygon(500, 150, 9));
  polygons.push(new Polygon(100, 300, 7));
  polygons.push(new Polygon(350, 400, 4));
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

  push()
  fill(200, 0, 200, 150);
  circle(mouseX, mouseY, 2 * mouseInteractionRadius);
  pop()

  dt = deltaTime
  
  
  // apply all external forces on all particles
  softBody.applyExternalForce(0, GRAVITY_FORCE_Y, dt);
  // softBody.applyExternalForce(10, 15, dt)
  softBody.applyExternalForce(mvx, mvy, dt);
  softBody.dampVelocity(DAMPING);
  
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
  // softBody.show();
  
  // if (collisionConstraints.length > 0) {
  //   console.log(collisionConstraints);
  //   for (let c of collisionConstraints) {
  //     let mass = softBody.masses[c.index]
  //     console.log(mass.position.x, mass.position.y);
  //     console.log(mass.predictedPosition.x, mass.predictedPosition.y);
  //     push()
  //     fill('blue')
  //     circle(mass.position.x, mass.position.y, 5)
  //     fill('red')
  //     circle(mass.predictedPosition.x, mass.predictedPosition.y, 5)
  //     pop()
  //   }
  //   noLoop()
  // }

  // let flag = false;
  // for (let mass of softBody.masses) {
  //   if (floor.collisionAlong(mass.position, mass.predictedPosition).collision) {
  //     flag = true;
  //   }
  // }
  // if (flag) {
  //   noLoop()
  // }

  //post-solve step
  softBody.updateVelocity(dt);
  softBody.updatePosition();

  softBody.updateCollidingMassVelocity(RESTITUTION, FRICTION);
}
