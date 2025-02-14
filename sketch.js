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
  softBody = new Soft_body(400, 300, 9, 5);
  floor = new Polygon([createVector(0, height), createVector(width, height), createVector(width, height - 20), createVector(0, height - 20)]);
  left = new Polygon([createVector(0, 0), createVector(20, 0), createVector(20, height), createVector(0, height)]);
  right = new Polygon([createVector(width - 20, 0), createVector(width, 0), createVector(width, height), createVector(width - 20, height)]);
  polygons.push(new Polygon(400, 400, 7));
  polygons.push(new Polygon(500, 250, 4));
  polygons.push(floor);
  polygons.push(left);
  polygons.push(right);

  s = new Spring(new Mass_point(100, 100, 1), new Mass_point(110, 250, 1), 30);
}

function draw() {
  background(51);
  // s.A.set_force_to_zero();
  // s.B.set_force_to_zero();
  // let gravity = s.A.calculate_gravity();
  // s.A.add_force(gravity);
  // s.B.add_force(gravity);
  // s.apply_spring_damp_force();
  // s.update_mass_point();
  // s.show();

  softBody.show();
  softBody.apply_forces();
  
  for (let polygon of polygons){
    polygon.show();
    softBody.handle_collision_with(polygon);
  }
  softBody.update();
  // noLoop()
}
