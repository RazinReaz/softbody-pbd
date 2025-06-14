let softBody;
let hashGrid;
let floor;
let left;
let right;
let polygons = [];

let pause = false;
let mvx, mvy;
let mouseInteractionRadius = 0;





function instructions() {
  push();
    noStroke();
    fill(201);
    textAlign(CENTER, TOP);
    textFont('Product Sans');
    textSize(15);    
    text("Click and drag with mouse to interact with the softbody", width / 2, 30);
  pop();
}

// user interaction
let showMethodButtons = [];
let rigMethodButtons = [];
let showMethods = [
  { name: "Rig", value: SoftBodyShowMethod.RIG },
  { name: "Surface", value: SoftBodyShowMethod.SURFACE }
];
let rigMethods = [
  { name: "Perimeter", value: SoftBodyRigMethod.PERIMETER },
  { name: "Grid", value: SoftBodyRigMethod.GRID }
];

let selectedRigMethod = SoftBodyRigMethod.GRID;
let selectedShowMethod = SoftBodyShowMethod.SURFACE;

function createSoftBodyUI() {
  // Show method buttons
  let showDiv = createDiv('Show Method:');
  showDiv.addClass('centered-div');
  showDiv.style('top', '10px');
  showMethodButtons = [];
  showMethods.forEach((method, idx) => {
    let btn = createButton(method.name);
    btn.parent(showDiv);
    btn.mousePressed(() => {
      selectedShowMethod = method.value;
      showMethodButtons.forEach(b => b.removeClass('selected'));
      btn.addClass('selected');
    });
    showMethodButtons.push(btn);
    if (idx === 0) btn.addClass('selected');
  });

  // Rig method buttons
  let rigDiv = createDiv('Rig Method:');
  rigDiv.addClass('centered-div');
  rigDiv.style('top', '60px');
  rigMethodButtons = [];
  rigMethods.forEach((method, idx) => {
    let btn = createButton(method.name);
    btn.parent(rigDiv);
    btn.mousePressed(() => {
      selectedRigMethod = method.value;
      rigMethodButtons.forEach(b => b.removeClass('selected'));
      btn.addClass('selected');
      restartSimulation();
    });
    rigMethodButtons.push(btn);
    if (idx === 0) btn.addClass('selected');
  });

  // Restart button
  let restartBtn = createButton('Restart Simulation');
  restartBtn.addClass('centered-div');
  restartBtn.style('top', '120px');
  restartBtn.mousePressed(() => {
    restartSimulation();
  });
}

function restartSimulation() {
  softBody = new Soft_body(350, 10, SOFTBODY_WIDTH, SOFTBODY_HEIGHT, SPRING_CONSTRAINT_STIFFNESS, SOLVER_ITERATIONS, selectedRigMethod, selectedShowMethod);
  hashGrid = new Grid(2 * SOFTBODY_SIZE, softBody.masses.length);
  // If you need to reset other state (like hashGrid), do it here as well
}

function setup() {
  createCanvas(800, 600);
  // UI stuff
  createSoftBodyUI();

  softBody = new Soft_body(350, 10, SOFTBODY_WIDTH, SOFTBODY_HEIGHT, SPRING_CONSTRAINT_STIFFNESS, SOLVER_ITERATIONS, selectedRigMethod, selectedShowMethod);
  hashGrid = new Grid(2 * SOFTBODY_SIZE, softBody.masses.length);

  floor = new Polygon([createVector(0, height - 50), createVector(width, height - 50), createVector(width, height - 50 - 50), createVector(0, height - 50 - 50)]);
  left = new Polygon([createVector(0, 0), createVector(50, height), createVector(100, height), createVector(50, 0)]);
  right = new Polygon([createVector(width, 0), createVector(width - 50, 0), createVector(width - 50, height), createVector(width, height)]);
  ceiling = new Polygon([createVector(0, 0), createVector(width, 0), createVector(width, -50), createVector(0, -50)]);
  polygons.push(new Polygon(500, 250, 9));
  // polygons.push(new Polygon(100, 300, 7));
  polygons.push(new Polygon(350, 400, 5));
  polygons.push(ceiling);
  polygons.push(floor);
  polygons.push(left);
  polygons.push(right);

}

let dt = 10;

function draw() {
  background(51);
  instructions();
  softBody.showMethod = selectedShowMethod;

  for (let mass of softBody.masses)
    hashGrid.update(mass)

  if (mouseIsPressed && mouseButton === LEFT) {
    mvx = (mouseX - pmouseX) * MOUSE_FORCE_MULTPLIER;
    mvy = (mouseY - pmouseY) * MOUSE_FORCE_MULTPLIER;
  } else {
    mvx = 0, mvy = 0;
  }
  
  // apply all external forces on all particles
  softBody.applyExternalForce(0, GRAVITY_FORCE_Y, dt);
  softBody.applyExternalForce(mvx, mvy, dt);
  softBody.dampVelocity(DAMPING);
  softBody.projectPosition(dt);
  
  for (let polygon of polygons){
    polygon.show();
  }
  softBody.show();

  softBody.generateCollisionConstraints(polygons);
  softBody.generateSelfCollisionConstraints(hashGrid);
  
  // solve step
  for (let i=0; i < SOLVER_ITERATIONS; i++) {
    // solve all constraints
    softBody.solveCollisionConstraints();
    softBody.solveSelfCollisionConstraints();
    softBody.solveSpringConstraints();
  }
  
  //post-solve step
  softBody.updateVelocity(dt);
  softBody.updateCollidingMassVelocity(RESTITUTION, FRICTION, dt);
  softBody.updatePosition();

}
