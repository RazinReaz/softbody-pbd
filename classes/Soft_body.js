class Soft_body {
  // for now the soft body is a recatangle.
  // we will build other shapes later
  constructor(positionx, positiony, width, height) {
    // create the structure
    this.springs = [];
    this.masses = [];
    let currentPos = {
      x : positionx,
      y : positiony
    }
    for (let i = 0; i < width; i++) {
      this.masses.push(new Mass_point(currentPos.x, currentPos.y, MASSPOINT_MASS));
      currentPos.x += SOFT_BODY_SIZE;
    } 
    currentPos.x -= SOFT_BODY_SIZE;
    
    for (let i = 0; i < height - 1; i++) {
      currentPos.y += SOFT_BODY_SIZE;
      this.masses.push(new Mass_point(currentPos.x, currentPos.y, MASSPOINT_MASS));
    } 
    
    for (let i = 0; i < width - 1; i++) {
      currentPos.x -= SOFT_BODY_SIZE;
      this.masses.push(new Mass_point(currentPos.x, currentPos.y, MASSPOINT_MASS));
    } 
    
    for (let i = 0; i < height - 2; i++) {
      currentPos.y -= SOFT_BODY_SIZE;
      this.masses.push(new Mass_point(currentPos.x, currentPos.y, MASSPOINT_MASS));
    } 

    // perimeter springs
    for (let i = 0; i < this.masses.length; i++) {
      let m1 = this.masses[i];
      let m2 = this.masses[(i+1) % this.masses.length];
      let d = dist(m1.position.x, m1.position.y, m2.position.x, m2.position.y);
      this.springs.push(new Spring(m1, m2, d))
    }

    for (let offset = width, k = 0; k < 2; offset++, k++){
      for (let i = 0; i < this.masses.length; i++) {
        let m1 = this.masses[i];
        let m2 = this.masses[(i+offset) % this.masses.length];
        let d = dist(m1.position.x, m1.position.y, m2.position.x, m2.position.y);
        this.springs.push(new Spring(m1, m2, d, false))
      }
    } 

    // rig springs

    // const top_left_index = 0;
    // const top_right_index = width - 1;
    // const bottom_right_index = top_right_index + height - 1;
    // const bottom_left_index = bottom_right_index + width - 1;
    
    // for (let i = top_left_index; i <= top_right_index - 1; i++) {
    //   let m1 = this.masses[i];
    //   let m2 = this.masses[bottom_left_index - 1 - i];
    //   let d = dist(m1.position.x, m1.position.y, m2.position.x, m2.position.y);
    //   this.springs.push(new Spring(m1, m2, d, false))
    // }
    // for (let i = top_left_index + 1; i <= top_right_index - 1; i++) {
    //   let m1 = this.masses[i];
    //   let m2 = this.masses[bottom_left_index - i];
    //   let d = dist(m1.position.x, m1.position.y, m2.position.x, m2.position.y);
    //   this.springs.push(new Spring(m1, m2, d, false))
    // }
    // for (let i = top_left_index + 1; i <= top_right_index; i++) {
    //   let m1 = this.masses[i];
    //   let m2 = this.masses[bottom_left_index - i + 1];
    //   let d = dist(m1.position.x, m1.position.y, m2.position.x, m2.position.y);
    //   this.springs.push(new Spring(m1, m2, d, false))
    // }

    // for (let i = bottom_left_index, k = 0; i <= bottom_left_index + height - 2; i++, k++) {
    //   let m1 = this.masses[i];
    //   let m2 = this.masses[bottom_right_index - k + 1];
    //   let d = dist(m1.position.x, m1.position.y, m2.position.x, m2.position.y);
    //   this.springs.push(new Spring(m1, m2, d, false))
    // }
    // let m1 = this.masses[0];
    // let m2 = this.masses[width];
    // let d = dist(m1.position.x, m1.position.y, m2.position.x, m2.position.y);
    // this.springs.push(new Spring(m1, m2, d, false))
    
    // for (let i = top_right_index, k = 0; i <= bottom_right_index - 1; i++, k++) {
    //   let m1 = this.masses[i];
    //   let m2 = this.masses[bottom_left_index + height - 2 - k];
    //   let d = dist(m1.position.x, m1.position.y, m2.position.x, m2.position.y);
    //   this.springs.push(new Spring(m1, m2, d, false))
    // }
    // for (let i = top_right_index + 1, k = 0; i <= bottom_right_index - 1; i++, k++) {
    //   let m1 = this.masses[i];
    //   let m2 = this.masses[bottom_left_index + height - 2 - k];
    //   let d = dist(m1.position.x, m1.position.y, m2.position.x, m2.position.y);
    //   this.springs.push(new Spring(m1, m2, d, false))
    // }

  }

  show() {
    for (let spring of this.springs) 
      spring.show();
    for (let mass of this.masses){
      mass.show();
      push();
      noFill();
      stroke(255);
      circle(mass.predictedPosition.x, mass.predictedPosition.y, 5);
      pop();
    } 
  }


  bounce(mass_point, reflection_point, normal) {
    // move the mass point
    mass_point.position = reflection_point;
    // update the mass point velocity
    let velocity = mass_point.velocity;
    let dot_product = velocity.dot(normal);
    let reflection_velocity = velocity.sub(normal.mult(2 * dot_product));
    mass_point.velocity = reflection_velocity; //! not sure
  }


  generatePolygonConstraints(polygons) {
    let collisionConstraints = [];
    for (let polygon of polygons) {
      for (let i = 0; i < this.masses.length; i++) {
        let mass = this.masses[i];
        let {collision, collidingLine, contact} = polygon.collisionAlong(mass.position, mass.predictedPosition)
        if (!collision)
          continue;
        // push();
        // fill('yellow');
        // circle(mass.position.x, mass.position.y, 2*mass.radius);
        // noFill();
        // stroke('blue')
        // circle(mass.predictedPosition.x, mass.predictedPosition.y, 3);
        // stroke(255)
        // circle(contact.x, contact.y, 3);
        // pop()
        // noLoop()

        let contactPoint = createVector(contact.x, contact.y)
        let n_hat = collidingLine.normal();
        
        collisionConstraints.push(
          {
            contactPoint : contactPoint,
            normal : n_hat,
            index: i,
          }
        )
      }
    }
    return collisionConstraints;
  }

  applyGravity(gravityConstant, deltaTime) {
    let gravityVector = createVector(0, gravityConstant)
    for (let mass of this.masses) {
      mass.velocity.add(p5.Vector.mult(gravityVector, deltaTime * mass.w));
    }
  }

  dampVelocity(dampingConstant) {
    for (let mass of this.masses) {
      mass.velocity.mult(dampingConstant);
    }
  }

  projectPosition(deltaTime) {
    // this function will calculate the predicted position of each of the particles
    for (let mass of this.masses) {
      mass.predictedPosition = p5.Vector.add(mass.position, p5.Vector.mult(mass.velocity, deltaTime));
    }
  }

  solveSpringConstraints() {
    for (let spring of this.springs) {
      let mass1 = spring.A;
      let mass2 = spring.B;
      let p1 = mass1.predictedPosition;
      let p2 = mass2.predictedPosition;

      let l = spring.resting_length;
      let vect = p5.Vector.sub(p2, p1);
      let d = vect.mag();

      let error = d - l;

      let v = p5.Vector.normalize(vect) // unit vector from p1 to p2 

      let lambda = -0.2 * error / (mass1.w + mass2.w);
      
      let delx1 = p5.Vector.mult(v, -lambda * mass1.w);
      let delx2 = p5.Vector.mult(v,  lambda * mass2.w);
      mass1.predictedPosition.add(delx1);
      mass2.predictedPosition.add(delx2)
    }
  }

  solveCollisionConstraints(constraints) {

    for (let constraint of constraints) {
      
      let {contactPoint, normal, index} = constraint;
      // evaluate the constraint C(p) = (p - con).n - d >= 0
      let mass = this.masses[index];
      let constraint_val = p5.Vector.sub(mass.predictedPosition, contactPoint).dot(normal) - mass.radius;

      if (constraint_val >= 0) continue;
      let deltax = p5.Vector.mult(normal, -1 * (constraint_val + EPSILON) * mass.w);
      mass.predictedPosition.add(deltax);


      // let velocity = mass.velocity;
      // let dot_product = velocity.dot(normal);
      // let reflection_velocity = velocity.sub(normal.mult(2 * dot_product));
      // mass.velocity = reflection_velocity.mult(BOUNCE_CONSTANT)

      // let dotProduct = mass.velocity.dot(normal);
      // let v_normal = p5.Vector.mult(normal, dotProduct);
      // let v_tangent = p5.Vector.sub(mass.velocity, v_normal)

      // mass.velocity = p5.Vector.mult(v_normal, BOUNCE_CONSTANT);
      // mass.velocity.add(p5.Vector.mult(v_tangent, SLIDE_CONSTANT));
    }
  }

  updateVelocity(deltaTIme) {
    for (let mass of this.masses) {
      let deltax = p5.Vector.sub(mass.predictedPosition, mass.position);
      mass.velocity = p5.Vector.div(deltax, deltaTIme);
    }
  }

  updatePosition() {
    for (let mass of this.masses) {
      mass.position = mass.predictedPosition;
    }
  }

  updateCollidingMassVelocity(constraints, restitution, friction) {
    for (let constraint of constraints) {
      let {contactPoint, normal, index} = constraint;
      // evaluate the constraint C(p) = (p - con).n - d >= 0
      let mass = this.masses[index];

      let dotProduct = mass.velocity.dot(normal);
      let v_normal = p5.Vector.mult(normal, dotProduct);
      let v_tangent = p5.Vector.sub(mass.velocity, v_normal);

      if (dotProduct < 0) {
        v_normal.mult(-restitution);
      }

      let len_v_t = v_tangent.mag()
      let frictionMagnitude = Math.max(0, 1 - friction * Math.abs(dotProduct) / len_v_t);
      v_tangent.mult(frictionMagnitude)

      mass.velocity = p5.Vector.add(v_normal, v_tangent);

      // mass.velocity = p5.Vector.mult(v_normal, BOUNCE_CONSTANT);
      // mass.velocity.add(p5.Vector.mult(v_tangent, SLIDE_CONSTANT));

      
    }
  }
}
