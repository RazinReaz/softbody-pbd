let softBody;
let floor;
let left;
let right;
let polygons = [];

let pause = false;

function setup() {
  createCanvas(800, 600);
  softBody = new Soft_body(450, 10, 15, 12);
  floor = new Polygon([createVector(0, height - 50), createVector(width, height - 50), createVector(width, height - 70), createVector(0, height - 70)]);
  left = new Polygon([createVector(0, 0), createVector(80, height), createVector(100, height), createVector(20, 0)]);
  right = new Polygon([createVector(width - 20, 0), createVector(width - 20, height), createVector(width, height), createVector(width, 0)]);
  polygons.push(new Polygon(500, 150, 9));
  polygons.push(new Polygon(200, 300, 7));
  polygons.push(new Polygon(500, 400, 5));
  polygons.push(floor);
  polygons.push(left);
  polygons.push(right);

}

function draw() {
  // frameRate(1);
  background(51);
  
  
  // apply all external forces on all particles
  softBody.applyGravity(GRAVITY_CONSTANT, deltaTime);
  softBody.dampVelocity(DAMPING);
  
  // calculate the predicted postions of all particles
  softBody.projectPosition(deltaTime);
  for (let polygon of polygons){
    polygon.show();
  }
  
  collisionConstraints = softBody.generatePolygonConstraints(polygons);
  // solve step
  for (let i=0; i < SOLVER_ITERATIONS; i++) {
    // solve all constraints
    softBody.solveSpringConstraints();
    softBody.solveCollisionConstraints(collisionConstraints);
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

  softBody.updateCollidingMassVelocity(collisionConstraints, RESTITUTION, FRICTION);
}
