let massIDCounter = 0;

class Mass_point{
    constructor(x, y, m){
        this.position = createVector(x, y);
        this.velocity = createVector(0, 0);
        this.w = 1/m;

        this.radius = MASSPOINT_RADIUS;
        this.predictedPosition = createVector(0,0);

        // for spatial hashing and self collision
        this.cellKey = -1;
        this.id = massIDCounter++;
        this.accumulatedVelocity = createVector(0, 0);
    }

    show(){
        push();
        // stroke(255);
        noStroke()
        fill(200,50,50);
        circle(this.position.x, this.position.y, 2 * this.radius);
        // circle(this.predictedPosition.x, this.predictedPosition.y, this.radius);
        pop();
    }

    generateKeyWith(otherMass){
        return this.id < otherMass.id ? `${this.id},${otherMass.id}` : `${otherMass.id},${this.id}`
    }

}