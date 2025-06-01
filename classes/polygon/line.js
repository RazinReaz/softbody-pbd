class Line{
    constructor(x1, y1, x2, y2){
        this.p = createVector(x1, y1);
        this.q = createVector(x2, y2);
        this.dist = dist(this.p.x, this.p.y, this.q.x, this.q.y);
        this.max_y = max(this.p.y, this.q.y);
        this.min_y = min(this.p.y, this.q.y);
        this.slope = (y2 - y1) / (x2 - x1);
        this.intercept = y1 - this.slope * x1;
    }

    show(strokeval = 255, weight = 1){
        push();
        strokeWeight(weight)
        stroke(strokeval);
        line(this.p.x, this.p.y, this.q.x, this.q.y);
        pop();
    }

    distance_from_point(point){
        let a = (this.q.x - this.p.x) * (this.p.y - point.y);
        let b = (this.p.x - point.x) * (this.q.y - this.p.y);
        return abs(a - b) / this.dist;
    }

    normal(){
        let edge = p5.Vector.sub(this.q, this.p);
        let normal = createVector(-edge.y, edge.x)
        return normal.normalize();
    }

    closest_point_on_line_from(point){
        let m = this.slope;
        let c = this.intercept;
        if ( m == 0 ) 
            return {
                x: point.x,
                y: this.p.y
            };

        if ( m == Infinity || m == -Infinity) 
            return {
                x: this.p.x,
                y: point.y
            }
        
        let b = point.y + point.x / m;
        let xp = (b - c) / (m + 1 / m);
        let yp = m * xp + c;
        return {
            x: xp,
            y: yp
        };
    }

    drawNormal(){
        // Midpoint of the line
        let mx = (this.p.x + this.q.x) / 2;
        let my = (this.p.y + this.q.y) / 2;
        // Get the normal vector (assume line_.normal() returns a p5.Vector)
        let n = this.normal();
        // Normalize and scale for visibility
        let n_scaled = n.copy().setMag(30);
        // Draw arrow for normal
        let x2 = mx + n_scaled.x;
        let y2 = my + n_scaled.y;
        // Draw the main line of the arrow
        line(mx, my, x2, y2);
        // Draw arrowhead
        push();
        translate(x2, y2);
        rotate(n.heading());
        line(0, 0, -7, -4);
        line(0, 0, -7, 4);
        pop();
    }
}