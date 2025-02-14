class Mass_point{
    constructor(x, y, m){
        this.position = createVector(x, y);
        this.velocity = createVector(0, 0);
        this.mass = m;
        this.force = createVector(0, 0);
        this.radius = MASSPOINT_RADIUS;
    }

    add_force(f){
        this.force.add(f);
    }

    calculate_gravity(){
        return createVector(0, GRAVITY_CONSTANT * this.mass);
    }

    set_force_to_zero(){
        this.force = createVector(0, 0);
    }

    show(){
        push();
        stroke(255);
        fill(200,50,50);
        circle(this.position.x, this.position.y, 5);
        pop();
    }

    update() {
        let acceleration = this.force.div(this.mass);
        let del_vel = p5.Vector.mult(acceleration, DELTA_TIME);
        this.velocity.add(del_vel);
        this.position.add(this.velocity);
    }
}