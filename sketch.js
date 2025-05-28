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
  softBody = new Soft_body(350, 10, 15, 13, SPRINT_CONSTRAINT_STIFFNESS, SOLVER_ITERATIONS);
  floor = new Polygon([createVector(0, height - 50), createVector(width, height - 50), createVector(width, height - 70), createVector(0, height - 70)]);
  left = new Polygon([createVector(0, 0), createVector(80, height), createVector(100, height), createVector(20, 0)]);
  right = new Polygon([createVector(width - 20, 0), createVector(width - 20, height), createVector(width, height), createVector(width, 0)]);
  polygons.push(new Polygon(500, 150, 9));
  polygons.push(new Polygon(100, 300, 7));
  polygons.push(new Polygon(500, 400, 5));
  polygons.push(floor);
  polygons.push(left);
  polygons.push(right);

}

function draw() {
  // frameRate(1);
  background(51);
  softBody.reset();
  if (mouseIsPressed && mouseButton === LEFT) {
    mouseInteractionRadius = MAX_MOUSE_INTERACTION_RADIUS;
    mvx = (mouseX - pmouseX) * MOUSE_FORCE_MULTPLIER;
    mvy = (mouseY - pmouseY) * MOUSE_FORCE_MULTPLIER;
  } else {
    mouseInteractionRadius = 0;
    mvx = 0;
    mvy = 0;
  }

  push()
  fill(200, 0, 200, 150);
  circle(mouseX, mouseY, 2 * mouseInteractionRadius);
  pop()
  
  
  // apply all external forces on all particles
  softBody.applyExternalForce(0, GRAVITY_FORCE_Y, deltaTime);
  // softBody.applyExternalForce(mvx, mvy, deltaTime)
  // softBody.applyMouseInteractionForce(mouseX, mouseY, mouseInteractionRadius, mvx, mvy, deltaTime)

  softBody.dampVelocity(DAMPING);
  
  // calculate the predicted postions of all particles
  softBody.projectPosition(deltaTime);
  for (let polygon of polygons){
    polygon.show();
  }
  
  softBody.generateCollisionConstraints(polygons);
  // solve step
  for (let i=0; i < SOLVER_ITERATIONS; i++) {
    // solve all constraints
    softBody.solveSpringConstraints();
    softBody.solveCollisionConstraints();
  }
  
  softBody.show();
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
  softBody.updateVelocity(deltaTime);
  softBody.updatePosition();

  softBody.updateCollidingMassVelocity(RESTITUTION, FRICTION);
}
