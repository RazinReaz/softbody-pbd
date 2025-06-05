
const SoftBodyShowMethod = {
  RIG: 0,
  SURFACE: 1,
  CONCAVE_HULL: 2,
  MARCHING_SQUARES: 3
  // Add more methods as needed in the future
};

const SoftBodyRigMethod = {
  DOUBLE: 0,
  PERIMETER: 1,
  GRID: 2
}


class Soft_body {
  // for now the soft body is a recatangle.
  // we will build other shapes later
  constructor(positionx, positiony, width, height, springStiffness, solverIterations) {
    this.springs = [];
    this.masses = [];
    this.width = width;
    this.height = height;
    this.collisionConstraints = [];
    this.selfCollisionConstraints = new Map(); // Map of [mass_point, mass_point] -> constraint data
    this.stiffnessMultiplier = 1 - Math.pow((1 - springStiffness), 1 / solverIterations)
    
    
    this.debug = false;
    this.showMethod = SoftBodyShowMethod.SURFACE;
    this.rigMethod = SoftBodyRigMethod.GRID;
    
    // for performance issues
    this._tmpContactPoint = createVector(0, 0);
    this._tempCollisionSubtract = createVector(0, 0);
    this._tempCollisionDeltaX = createVector(0, 0);
    this._tempSelfCollisionTotalDeltaV = createVector(0, 0);
    this._forceVector = createVector(0, 0);
    this._tempInGeneral = createVector(0, 0);
    this._tempInGeneral2 = createVector(0, 0);
    
    let currentPos = {
      x : positionx,
      y : positiony
    }

    if (this.rigMethod === SoftBodyRigMethod.PERIMETER) {

      for (let i = 0; i < width; i++) {
        this.masses.push(new Mass_point(currentPos.x, currentPos.y, MASSPOINT_MASS));
        currentPos.x += SOFTBODY_SIZE;
      } 
      currentPos.x -= SOFTBODY_SIZE;
      
      
      for (let i = 0; i < height - 1; i++) {
        currentPos.y += SOFTBODY_SIZE;
        this.masses.push(new Mass_point(currentPos.x, currentPos.y, MASSPOINT_MASS));
      } 
      
      for (let i = 0; i < width - 1; i++) {
        currentPos.x -= SOFTBODY_SIZE;
        this.masses.push(new Mass_point(currentPos.x, currentPos.y, MASSPOINT_MASS));
      } 
      
      for (let i = 0; i < height - 2; i++) {
        currentPos.y -= SOFTBODY_SIZE;
        this.masses.push(new Mass_point(currentPos.x, currentPos.y, MASSPOINT_MASS));
      } 
  
      // perimeter springs
      for (let i = 0; i < this.masses.length; i++) {
        let m1 = this.masses[i];
        let m2 = this.masses[(i+1) % this.masses.length];
        let d = dist(m1.position.x, m1.position.y, m2.position.x, m2.position.y);
        this.springs.push(new Spring(m1, m2, d, this.stiffnessMultiplier))
      }
  
      for (let offset = width, k = 0; k < 2; offset++, k++){
        for (let i = 0; i < this.masses.length; i++) {
          let m1 = this.masses[i];
          let m2 = this.masses[(i+offset+k) % this.masses.length];
          let d = dist(m1.position.x, m1.position.y, m2.position.x, m2.position.y);
          this.springs.push(new Spring(m1, m2, d, this.stiffnessMultiplier * (width * SOFTBODY_SIZE / d) * 0.3)) //! experimental
        }
      } 

    } else if (this.rigMethod === SoftBodyRigMethod.DOUBLE) {
      for (let i = 0; i < width; i++) {
        this.masses.push(new Mass_point(currentPos.x, currentPos.y, MASSPOINT_MASS));
        currentPos.x += SOFTBODY_SIZE;
      } 
      currentPos.x -= SOFTBODY_SIZE;
      for (let i = 0; i < height - 1; i++) {
        currentPos.y += SOFTBODY_SIZE;
        this.masses.push(new Mass_point(currentPos.x, currentPos.y, MASSPOINT_MASS));
      } 
      for (let i = 0; i < width - 1; i++) {
        currentPos.x -= SOFTBODY_SIZE;
        this.masses.push(new Mass_point(currentPos.x, currentPos.y, MASSPOINT_MASS));
      } 
      for (let i = 0; i < height - 2; i++) {
        currentPos.y -= SOFTBODY_SIZE;
        this.masses.push(new Mass_point(currentPos.x, currentPos.y, MASSPOINT_MASS));
      } 
      
      // perimeter springs
      for (let i = 0; i < this.masses.length; i++) {
        let m1 = this.masses[i];
        let m2 = this.masses[(i+1) % this.masses.length];
        let d = dist(m1.position.x, m1.position.y, m2.position.x, m2.position.y);
        this.springs.push(new Spring(m1, m2, d, this.stiffnessMultiplier))
      }
      // perimeter springs strengthening
      for (let offset = 5, k = 0; k < 3; offset++, k++){
        for (let i = 0; i < this.masses.length; i++) {
          let m1 = this.masses[i];
          let m2 = this.masses[(i+offset + k) % this.masses.length];
          let d = dist(m1.position.x, m1.position.y, m2.position.x, m2.position.y);
          this.springs.push(new Spring(m1, m2, d, this.stiffnessMultiplier))
        }
      } 
      let innerSize = 0.5;
      const innerBodySize = SOFTBODY_SIZE * innerSize
      currentPos = {
        x: positionx + width * (1 - innerSize) * 0.5 * SOFTBODY_SIZE,
        y: positiony + height * (1 - innerSize) * 0.5 * SOFTBODY_SIZE,
      }
      let tempOffset = this.masses.length;
      for (let i = 0; i < width; i++) {
        this.masses.push(new Mass_point(currentPos.x, currentPos.y, MASSPOINT_MASS));
        currentPos.x += innerBodySize;
      } 
      currentPos.x -= innerBodySize;
      
      for (let i = 0; i < height - 1; i++) {
        currentPos.y += innerBodySize;
        this.masses.push(new Mass_point(currentPos.x, currentPos.y, MASSPOINT_MASS));
      } 
      
      for (let i = 0; i < width - 1; i++) {
        currentPos.x -= innerBodySize;
        this.masses.push(new Mass_point(currentPos.x, currentPos.y, MASSPOINT_MASS));
      } 
      
      for (let i = 0; i < height - 2; i++) {
        currentPos.y -= innerBodySize;
        this.masses.push(new Mass_point(currentPos.x, currentPos.y, MASSPOINT_MASS));
      } 
      // inner springs
      for (let i = tempOffset; i < this.masses.length; i++) {
        let m1 = this.masses[i];
        let m2 = this.masses[(i+1) % this.masses.length];
        let d = dist(m1.position.x, m1.position.y, m2.position.x, m2.position.y);
        this.springs.push(new Spring(m1, m2, d, this.stiffnessMultiplier))
      }
      // inner springs strengthening
      for (let offset = 3, k = 0; k < 5; offset++, k++){
        for (let i = 0; i < this.masses.length - tempOffset; i++) {
          let m1 = this.masses[tempOffset + i];
          let m2 = this.masses[(tempOffset + i + offset + k) % tempOffset + tempOffset];
          let d = dist(m1.position.x, m1.position.y, m2.position.x, m2.position.y);
          this.springs.push(new Spring(m1, m2, d, this.stiffnessMultiplier))
        }
      } 
      // inner-outer springs 
      for (let i = tempOffset, j = 0; i < this.masses.length; i++, j++) {
        let inner = this.masses[i];
        let outer1 = this.masses[(j-2 + tempOffset)%tempOffset]; 
        let outer2 = this.masses[j]; 
        let outer3 = this.masses[(j+2)%tempOffset]; 
        let d1 = dist(inner.position.x, inner.position.y, outer1.position.x, outer1.position.y);
        let d2 = dist(inner.position.x, inner.position.y, outer2.position.x, outer2.position.y);
        let d3 = dist(inner.position.x, inner.position.y, outer3.position.x, outer3.position.y);
        this.springs.push(new Spring(inner, outer1, d1, this.stiffnessMultiplier))
        this.springs.push(new Spring(inner, outer2, d2, this.stiffnessMultiplier))
        this.springs.push(new Spring(inner, outer3, d3, this.stiffnessMultiplier))
      }
    } else if (this.rigMethod === SoftBodyRigMethod.GRID) {
      for (let i = 0; i < height; i++) {
        for (let j = 0; j < width; j++) {
          this.masses.push(new Mass_point(currentPos.x, currentPos.y, MASSPOINT_MASS));
          currentPos.x += SOFTBODY_SIZE;
        }
        currentPos.x = positionx;
        currentPos.y += SOFTBODY_SIZE;
      }

      for (let i = 0; i < height - 1; i++) {
        for (let j = 0; j < width - 1; j++) {
          let idx = i * width + j;
          let m1 = this.masses[idx];
          let m2 = this.masses[idx + 1];
          let m3 = this.masses[idx + width];
          let m4 = this.masses[idx + width + 1];
          let d12 = dist(m1.position.x, m1.position.y, m2.position.x, m2.position.y);
          let d13 = dist(m1.position.x, m1.position.y, m3.position.x, m3.position.y);
          let d14 = dist(m1.position.x, m1.position.y, m4.position.x, m4.position.y);
          let d23 = dist(m2.position.x, m2.position.y, m3.position.x, m3.position.y);
          this.springs.push(new Spring(m1, m2, d12, this.stiffnessMultiplier));
          this.springs.push(new Spring(m1, m3, d13, this.stiffnessMultiplier));
          this.springs.push(new Spring(m1, m4, d14, this.stiffnessMultiplier));
          this.springs.push(new Spring(m2, m3, d23, this.stiffnessMultiplier));

          if ( j === width - 2) {
            let d24 = dist(m2.position.x, m2.position.y, m4.position.x, m4.position.y);
            this.springs.push(new Spring(m2, m4, d24, this.stiffnessMultiplier));
          }
          if (i === height - 2) {
            let d34 = dist(m3.position.x, m3.position.y, m4.position.x, m4.position.y);
            this.springs.push(new Spring(m3, m4, d34, this.stiffnessMultiplier));
          }
        }
      } 


    }

  }

  #getPerimeterPoints() {
    if (this.rigMethod === SoftBodyRigMethod.PERIMETER) {
      return this.masses
    }

    else if (this.rigMethod === SoftBodyRigMethod.GRID) {
      let perimeterPoints = [];
      // push the first row into the perimeter
      for (let i = 0; i < this.width; i++) {
        perimeterPoints.push(this.masses[i]);
      }
      // push the last column into the perimeter
      for (let i = 2 * this.width - 1; i < this.masses.length; i += this.width) {
        perimeterPoints.push(this.masses[i]);
      }
      // push the last row into the perimeter in reverse order
      for (let i = this.masses.length - 1; i > this.masses.length - this.width; i--) {
        perimeterPoints.push(this.masses[i]);
      }
      // push the first column into the perimeter in reverse order
      for (let i = this.masses.length - this.width; i > 0; i -= this.width) {
        perimeterPoints.push(this.masses[i]);
      } 
      return perimeterPoints;
    }
  }

  show() {
    if (this.showMethod === SoftBodyShowMethod.RIG) {
      for (let spring of this.springs) 
        spring.show();
      for (let mass of this.masses)
        mass.show();
    }
    else if (this.showMethod === SoftBodyShowMethod.SURFACE) {
      // convex hull
      // let points = grahamScan(this.masses.map(m => m.position));
      const perimeterPoints = this.#getPerimeterPoints();
      let points = perimeterPoints.map(mass => mass.position);
      push();
      fill(200, 50, 50, 200);
      beginShape();
      for (let point of points) {
        vertex(point.x, point.y);
      }
      endShape(CLOSE);
      pop();
      // noLoop();
    }

    if (!this.debug) return;

    for (let constraint of this.collisionConstraints) {
      let mass = this.masses[constraint.index];
      let p = constraint.contactPoint
      push();
      fill('blue')
      circle(mass.position.x, mass.position.y, 2*mass.radius);
      circle(mass.predictedPosition.x, mass.predictedPosition.y, mass.radius);
      fill('white')
      circle(p.x, p.y, 3);
      pop();
    } 
  }


  generateCollisionConstraints(polygons) {
    this.collisionConstraints.length = 0;
    for (let i = 0; i < this.masses.length; i++) {
      let mass = this.masses[i];
      let minT = Infinity;
      let contact = null;
      let collidingLine = null;

      for (const polygon of polygons) {
        const {collision, collidingLine: line, contact: contactCandidate, t} = polygon.collisionAlongRay(mass.position, mass.predictedPosition)
        if (!collision) continue;

        if (this.debug) {
          console.log(t == -1 && collision ? "Static Collision" : "Continuous collision");
        }

        if (t < minT) {
          minT = t;
          contact = contactCandidate;
          collidingLine = line;
        }
      }
      if (minT === Infinity) continue;


      // collision found between mass and a polygon
      this._tmpContactPoint.set(contact.x, contact.y);
      let contactPoint = this._tmpContactPoint.copy(); //! copy or not?
      let normal = collidingLine.normal();
      
      this.collisionConstraints.push({
        contactPoint,
        normal,
        index: i,
      });
    }
  }

  generateSelfCollisionConstraints(hashGrid) {
    // Reset the map
    this.selfCollisionConstraints.clear();

    for (const mass of this.masses) {
        const neighbours = hashGrid.query(mass.predictedPosition);

        for (const neighbour of neighbours) {
            // Skip self
            // if (neighbor === mass) continue;

            const key = mass.generateKeyWith(neighbour)
            const distSq = distSquared(mass.predictedPosition, neighbour.predictedPosition);
            if (distSq == 0) continue; 
            const sum_radius = mass.radius + neighbour.radius;
            if (distSq < sum_radius * sum_radius && !this.selfCollisionConstraints.has(key)) {
                // Add the constraint to the map
                this.selfCollisionConstraints.set(key, {
                    massA: mass,
                    massB: neighbour,
                });
            }
        }
    }
  }

  applyExternalForce(extForceX, extForceY, deltaTime) {
    // let forceVector = createVector(extForceX, extForceY);
    for (let mass of this.masses) {
      this._forceVector.set(extForceX, extForceY);
      mass.velocity.add(this._forceVector.mult(deltaTime * mass.w));
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
      this._tempInGeneral.set(mass.velocity.x, mass.velocity.y);
      this._tempInGeneral.mult(deltaTime);
      this._tempInGeneral.add(mass.position);
      mass.predictedPosition.set(this._tempInGeneral.x, this._tempInGeneral.y);
      // mass.predictedPosition = p5.Vector.add(mass.position, p5.Vector.mult(mass.velocity, deltaTime));
    }
  }

  #solveSingleSpringConstraint(spring) {
    const mass1 = spring.A;
    const mass2 = spring.B;
    const p1 = mass1.predictedPosition;
    const p2 = mass2.predictedPosition;
    const l = spring.resting_length;
    const vect = p5.Vector.sub(p2, p1);
    const d = vect.mag();
    
    // Early exit if masses are coincident or both are immovable
    if (d === 0 || (mass1.w + mass2.w === 0)) return;
    
    const constraintValue = d - l;
    const n = vect.copy().normalize() // unit vector from p1 to p2 
    const lambda = -1.0 * constraintValue / (mass1.w + mass2.w);
    
    // const delx1 = p5.Vector.mult(n, -lambda * mass1.w * spring.stiffness);
    // const delx2 = p5.Vector.mult(n,  lambda * mass2.w * spring.stiffness);
    // mass1.predictedPosition.add(delx1);
    // mass2.predictedPosition.add(delx2);
    this._tempInGeneral.set(n.x, n.y).mult(-lambda * mass1.w * spring.stiffness);
    mass1.predictedPosition.add(this._tempInGeneral);
    this._tempInGeneral.set(n.x, n.y).mult(lambda * mass2.w * spring.stiffness);
    mass2.predictedPosition.add(this._tempInGeneral);
  }

  solveSpringConstraints() {
    for (let spring of this.springs) {
      this.#solveSingleSpringConstraint(spring);
    }
  }

  solveCollisionConstraints() {
    for (let constraint of this.collisionConstraints) {
      
      let {contactPoint, normal, index} = constraint;
      // evaluate the constraint C(p) = (p - con).n - d >= 0
      let mass = this.masses[index];
      this._tempCollisionSubtract.set(mass.predictedPosition.x, mass.predictedPosition.y).sub(contactPoint);
      let constraint_val = this._tempCollisionSubtract.dot(normal) - mass.radius;

      if (constraint_val >= 0) continue;
      this._tempCollisionDeltaX.set(normal.x, normal.y).mult(-1 * (constraint_val + EPSILON)* mass.w);
      mass.predictedPosition.add(this._tempCollisionDeltaX);
      // let deltax = p5.Vector.mult(normal, -1 * (constraint_val + EPSILON) * mass.w);
      // mass.predictedPosition.add(deltax);
    }
  }

  solveSelfCollisionConstraints(){
    // iterate through the map
    for (const [key, constraint] of this.selfCollisionConstraints.entries()) {
      let { massA, massB } = constraint;
      // n : A => B
      this._tempCollisionSubtract.set(massB.predictedPosition.x, massB.predictedPosition.y).sub(massA.predictedPosition);
      const distance = this._tempCollisionSubtract.mag();
      const constraintValue = distance - massA.radius - massB.radius;

      if (constraintValue >= 0) continue;

      const n = this._tempCollisionSubtract.normalize();
      const lambda = -1.0 * constraintValue / (massA.w + massB.w);
      const stiffness = 1;
      this._tempInGeneral.set(n.x, n.y).mult(-lambda * massA.w * stiffness);
      massA.predictedPosition.add(this._tempInGeneral);
      this._tempInGeneral.set(n.x, n.y).mult(lambda * massB.w * stiffness);
      massB.predictedPosition.add(this._tempInGeneral);
    }
  }

  updateVelocity(deltaTIme) {
    for (let mass of this.masses) {
      // let deltax = p5.Vector.sub(mass.predictedPosition, mass.position);
      this._tempInGeneral.set(mass.predictedPosition.x, mass.predictedPosition.y);
      this._tempInGeneral.sub(mass.position).div(deltaTIme);
      // mass.velocity = p5.Vector.div(deltax, deltaTIme);
      mass.velocity.set(this._tempInGeneral.x, this._tempInGeneral.y);
    }
  }

  updatePosition() {
    for (let mass of this.masses) {
      // mass.position = mass.predictedPosition.copy();
      mass.position.set(mass.predictedPosition.x, mass.predictedPosition.y);
    }
  }

  updateCollidingMassVelocity(restitution, friction) {
    for (let constraint of this.collisionConstraints) {
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
      // friction magnitude is some number between 0 and 1
      // the more the velocity is aligned to the normal, the more the friction magnitude
      v_tangent.mult(frictionMagnitude)

      mass.velocity = p5.Vector.add(v_normal, v_tangent);

      // mass.velocity = p5.Vector.mult(v_normal, BOUNCE_CONSTANT);
      // mass.velocity.add(p5.Vector.mult(v_tangent, SLIDE_CONSTANT));      
    }

    // do something for self collisions too
    // for (const [key, constraint] of this.selfCollisionConstraints.entries()) {
    //   let { massA, massB } = constraint;
    //   massA.velocity.set(0, 0);
    //   massB.velocity.set(0, 0);
    // }
    for (let mass of this.masses){
      mass.accumulatedVelocity.set(0, 0);
    }

    for (const [key, constraint] of this.selfCollisionConstraints.entries()) {
      let { massA:mass1, massB:mass2 } = constraint;
      let n = p5.Vector.sub(mass1.predictedPosition, mass2.predictedPosition);
      let relVel = p5.Vector.sub(mass1.velocity, mass2.velocity);
      const dotP = n.dot(relVel);
      const sumW = mass1.w + mass2.w;
      const dsq = n.magSq();
      const m = 2 * dotP / (sumW * dsq)

      this._tempInGeneral.set(n.x, n.y).mult(m * mass1.w * -1);
      mass1.accumulatedVelocity.add(this._tempInGeneral);
      this._tempInGeneral.set(n.x, n.y).mult(m * mass2.w);
      mass2.accumulatedVelocity.add(this._tempInGeneral);
    }

    for (let mass of this.masses){
      mass.velocity.add(mass.accumulatedVelocity);
    }
  }
}