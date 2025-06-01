class Spring{
    constructor(mass_pointA, mass_pointB, resting_length, stiffness, visible = true){
        this.A = mass_pointA;
        this.B = mass_pointB;
        this.resting_length = resting_length;
        this.stiffness = stiffness;
        this.visible = visible;
    }
    show(){
        if (!this.visible) return;

        push()
        stroke(200)
        line(this.A.position.x, this.A.position.y, this.B.position.x, this.B.position.y);
        pop()
    }

    get_length() {
        return p5.Vector.sub(this.A.position, this.B.position).mag()
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

}