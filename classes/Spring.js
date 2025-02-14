class Spring{
    constructor(mass_pointA, mass_pointB, resting_length, spring_constant = SPRING_SPRING_CONSTANT, damping_constant = SPRING_DAMPING_CONSTANT){
        this.A = mass_pointA;
        this.B = mass_pointB;
        this.resting_length = resting_length;
        this.spring_constant = spring_constant;
        this.damping_constant = damping_constant;
    }
    show(){
        // let n = 7;
        // let t = p5.Vector.sub(this.B.position, this.A.position).normalize();
        // let normal = createVector(-t.y, t.x);
        // push();
        //     noFill();
        //     stroke(255);
        //     beginShape()
        //     vertex(this.A.position.x, this.A.position.y);
        //     for(let i = 1; i < n - 1; i++){
        //         let m = pow(-1, i);
        //         let p = p5.Vector.lerp(this.A.position, this.B.position, i / n);
        //         let point = p5.Vector.add(p, p5.Vector.mult(normal, m * 10));
        //         vertex(point.x, point.y);
        //     }
        //     vertex(this.B.position.x, this.B.position.y);
        //     endShape();
        // pop();


        push()
        stroke(200)
        line(this.A.position.x, this.A.position.y, this.B.position.x, this.B.position.y);
        pop()
        this.A.show();
        this.B.show();
    }

    get_length() {
        return p5.Vector.sub(this.A.position, this.B.position).mag()
    }

    apply_spring_damp_force(){
        // spring force
        let spring_force_scalar = this.calculate_spring_force();
        let damping_force_scalar = this.calculate_damping_force();
        // console.log(spring_force_scalar, damping_force_scalar)
        let total_force_scalar = spring_force_scalar + damping_force_scalar;
        let direction = p5.Vector.sub(this.B.position, this.A.position).normalize();
        let A_force = p5.Vector.mult(direction, total_force_scalar)
        let B_force = p5.Vector.mult(direction, -total_force_scalar)
        this.A.add_force(A_force);
        this.B.add_force(B_force);
    }

    apply_self_collision() {
        let distance_vec = p5.Vector.sub(this.A.position, this.B.position);
        let distanceSq = distance_vec.magSq();
        let min_dist = this.resting_length * 0.2;
        if (distanceSq > min_dist * min_dist || distanceSq == 0)
            return;

        // the closer the distance, the more the repulsive force along the normal
        let distance = Math.sqrt(distanceSq);
        let normal = distance_vec.copy().div(distance);
        let scalar = (distance - min_dist/2) * 100;
        let repulsive_force = normal.mult(scalar);

        this.A.add_force(repulsive_force);
        this.B.add_force(repulsive_force.mult(-1))
    }


    update_mass_point(){
        /**
         * update both mass points according to the force applied to them
         */
        this.A.update();
        this.B.update();
    }

    calculate_spring_force(){
        /**
         * F = kx
         * @returns {number} spring force scalar
         */
        let distance = p5.Vector.sub(this.B.position, this.A.position).mag();
        return this.spring_constant * (distance - this.resting_length);
    }

    calculate_damping_force(){
        /**
         * F = -bv
         * @returns {number} damping force scalar
         */
        let relative_velocity = p5.Vector.sub(this.B.velocity, this.A.velocity);
        let relative_velocity_scalar = relative_velocity.dot(p5.Vector.sub(this.B.position, this.A.position).normalize());
        return this.damping_constant * relative_velocity_scalar;
    }


}