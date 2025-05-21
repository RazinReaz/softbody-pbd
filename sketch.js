let softBody;
let floor;
let left;
let right;
let polygons = [];
let s;

let pause = false;



function setup() {
  createCanvas(800, 600);
  // frameRate(5);
  softBody = new Soft_body(400, 100, 6, 5);
  floor = new Polygon([createVector(0, height), createVector(width, height), createVector(width, height - 20), createVector(0, height - 20)]);
  left = new Polygon([createVector(0, 0), createVector(20, 0), createVector(20, height), createVector(0, height)]);
  right = new Polygon([createVector(width - 20, 0), createVector(width, 0), createVector(width, height), createVector(width - 20, height)]);
  polygons.push(new Polygon(300, 400, 7));
  polygons.push(new Polygon(500, 250, 4));
  polygons.push(floor);
  polygons.push(left);
  polygons.push(right);

  // s = new Spring(new Mass_point(100, 100, 1), new Mass_point(110, 250, 1), 30);
}

function draw() {
  background(51);
  softBody.show();
  // noLoop();
  
  // apply all external forces on all particles
  softBody.applyGravity(GRAVITY_CONSTANT, deltaTime);
  
  // calculate the predicted postions of all particles
  softBody.projectPosition(deltaTime);

  collisionConstraints = softBody.generatePolygonConstraints(polygons);

  // solve step
  for (let i=0; i < SOLVER_ITERATIONS; i++) {
    // solve all constraints
    softBody.solveSpringConstraints();
    softBody.solveCollisionConstraints(collisionConstraints);
  }
  

  //post-solve step
  softBody.updateVelocity(deltaTime);
  softBody.updatePosition();
  for (let polygon of polygons){
    polygon.show();
  }
  // softBody.update();
  // noLoop()
}
