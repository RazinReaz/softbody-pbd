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
        if ( m == 0 ) return createVector(point.x, this.p.y);
        if ( m == Infinity || m == -Infinity) return createVector(this.p.x, point.y);
        
        let b = point.y + point.x / m;
        let xp = (b - c) / (m + 1 / m);
        let yp = m * xp + c;
        return createVector(xp, yp);
    }
}