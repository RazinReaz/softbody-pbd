class Soft_body {
  // for now the soft body is a recatangle.
  // we will build other shapes later
  constructor(positionx, positiony, width, height) {
    // create the structure
    this.springs = [];
    this.grid = [];
    for (let i = 0; i < height; i++) {
      this.grid.push([]);
      for (let j = 0; j < width; j++) {
        this.grid[i].push(
          new Mass_point(
            positionx + j * SOFT_BODY_SIZE,
            positiony + i * SOFT_BODY_SIZE,
            10
          )
        );
      }
    }
    for (let i = 0; i < height; i++) {
      for (let j = 0; j < width; j++) {
        if (i != height - 1)
          this.springs.push(
            new Spring(this.grid[i][j], this.grid[i + 1][j], SOFT_BODY_SIZE)
          );
        if (j != width - 1)
          this.springs.push(
            new Spring(this.grid[i][j], this.grid[i][j + 1], SOFT_BODY_SIZE)
          );
        if (i != height - 1 && j != width - 1)
          this.springs.push(
            new Spring(this.grid[i][j], this.grid[i + 1][j + 1], SOFT_BODY_SIZE * 1.414)
          );
        if (i != 0 && j != width - 1)
          this.springs.push(
            new Spring(this.grid[i][j], this.grid[i - 1][j + 1], SOFT_BODY_SIZE * 1.414)
          );
      }
    }
  }

  show() {
    for (let i = 0; i < this.springs.length; i++) {
      this.springs[i].show();
    }
  }

  apply_forces() {
    for (let i = 0; i < this.grid.length; i++) {
      for (let j = 0; j < this.grid[i].length; j++) {
        this.grid[i][j].set_force_to_zero();
        let grav = this.grid[i][j].calculate_gravity();
        this.grid[i][j].add_force(grav);
      }
    }

    for (let spring of this.springs) {
      spring.apply_spring_damp_force();
      spring.apply_self_collision();
    }
  }

  update() {
    for (let i = 0; i < this.grid.length; i++) {
      for (let j = 0; j < this.grid[i].length; j++) {
        let acceleration = this.grid[i][j].force.div(this.grid[i][j].mass);
        let del_vel = p5.Vector.mult(acceleration, DELTA_TIME)
        this.grid[i][j].velocity.add(del_vel);
        this.grid[i][j].position.add(this.grid[i][j].velocity);
      }
    }
  }

  bounce(mass_point, reflection_point, normal) {
    // move the mass point
    mass_point.position = reflection_point;
    // update the mass point velocity
    let velocity = mass_point.velocity;
    let dot_product = velocity.dot(normal);
    let reflection_velocity = velocity.sub(normal.mult(2 * dot_product));
    mass_point.velocity = reflection_velocity;
  }

  handle_collision_with(polygon) {
    for (let i = 0; i < this.grid.length; i++) {
      for (let j = 0; j < this.grid[i].length; j++) {
        if (polygon.collides_with(this.grid[i][j])){
          let { point, normal } = polygon.get_reflection_point_and_normal(this.grid[i][j]);
          this.bounce(this.grid[i][j], point, normal);
        }
      }
    }
  }


}
