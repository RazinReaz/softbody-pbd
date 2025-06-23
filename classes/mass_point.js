let massIDCounter = 0;

class Mass_point{
    constructor(x, y, m){
        this.position = createVector(x, y);
        this.velocity = createVector(0, 0);
        this.w = 1/m;

        this.radius = MASSPOINT_RADIUS;
        this.repulseRadius = this.radius * REPULSE_RADIUS_MULTIPLIER;
        this.predictedPosition = createVector(0,0);

        // for spatial hashing and self collision
        this.cellKey = -1;
        this.id = massIDCounter++;
        this.accumulatedVelocity = createVector(0, 0);
    }

    show(r=200, g=50, b=50){
        push();
        // stroke(255);
        noStroke()
        fill(r, g, b);
        circle(this.position.x, this.position.y, 2 * this.radius);
        // circle(this.predictedPosition.x, this.predictedPosition.y, this.radius);
        pop();
    }

    generateKeyWith(otherMass){
        return this.id < otherMass.id ? `${this.id},${otherMass.id}` : `${otherMass.id},${this.id}`
    }

}