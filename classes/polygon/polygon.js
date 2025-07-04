function polygon_points(x, y, n, radius){
    let points = [];
    let a = 2 * PI / n;
    for (let i = 0; i < n; i++){
        let r =  random(1,2);
        points.push(createVector(x + r * radius * cos(i * a), y - r * radius * sin(i * a)));
    }
    return points;
}

// Checks if point c is on segment ab
function onSegment(a, b, c) {
    return Math.min(a.x, b.x) <= c.x && c.x <= Math.max(a.x, b.x) && Math.min(a.y, b.y) <= c.y && c.y <= Math.max(a.y, b.y);
}

function getLineIntersection(p0, p1, p2, p3) {
    //! this function is not tested. it was copied from stack overflow
    let s10 = {
        x: p1.x - p0.x,
        y: p1.y - p0.y
    }
    let s32 = {
        x: p3.x - p2.x,
        y: p3.y - p2.y
    }

    let denom = s10.x * s32.y - s32.x * s10.y;
    if (denom == 0) return {collision: false, contactPoint: null, t : Infinity}; // colinear

    let denomPositive = denom > 0;

    let s02 = {
        x: p0.x - p2.x,
        y: p0.y - p2.y
    }

    let s_numer = s10.x * s02.y - s10.y * s02.x;
    if ((s_numer < 0) == denomPositive)
        return {collision: false, contactPoint: null, t : Infinity}; // no collision

    let t_numer = s32.x * s02.y - s32.y * s02.x;
    if ((t_numer < 0) === denomPositive)
        return { collision: false, contactPoint: null, t : Infinity }; // No collision

    if (((s_numer > denom) === denomPositive) || ((t_numer > denom) === denomPositive))
        return { collision: false, contactPoint: null, t : Infinity }; // No collision

    // Collision detected
    let t_ = t_numer / denom;
    let intersection = {
        x: p0.x + (t_ * s10.x),
        y: p0.y + (t_ * s10.y)
    };
    // console.log("intersection: ", intersection)
    // console.log("p0.x + (t * s10.x): ", p0.x + (t * s10.x))
    // console.log("p0.y + (t * s10.y): ", p0.y + (t * s10.y))

    return { collision: true, contactPoint: intersection, t: t_ };
}



class Polygon{
    constructor(x, y, n, radius = 50){
        this.points = [];
        if (Array.isArray(x)) 
        {
            // Sort points counter-clockwise
            // Compute centroid
            let cx = 0, cy = 0;
            for (let pt of x) {
                cx += pt.x;
                cy += pt.y;
            }
            cx /= x.length;
            cy /= x.length;
            // Sort by angle from centroid
            let sorted = x.slice().sort((a, b) => {
                let angleA = Math.atan2(a.y - cy, a.x - cx);
                let angleB = Math.atan2(b.y - cy, b.x - cx);
                return angleB - angleA;
            });
            for (let pt of sorted) {
                this.points.push(pt);
            }
        } 
        else 
        {
            this.points = polygon_points(x, y, n, radius);
        }
        this.lines = [];
        for (let i = 0; i < this.points.length; i++)
        {
            this.lines.push(new Line(this.points[i].x, this.points[i].y, this.points[(i+1)%this.points.length].x, this.points[(i+1)%this.points.length].y));
        }
        

        //bounding box
        this.max_x = 0;
        this.max_y = 0;
        this.min_x = 10000;
        this.min_y = 10000;
        for (let point of this.points){
            if (point.x > this.max_x) this.max_x = point.x;
            if (point.x < this.min_x) this.min_x = point.x;
            if (point.y > this.max_y) this.max_y = point.y;
            if (point.y < this.min_y) this.min_y = point.y;
        }
    }

    show(fillColor = 250){
        // TODO fill using other methods
        push()
        fill(fillColor)
        beginShape();
        for (let point of this.points){
            vertex(point.x, point.y);
        }
        endShape(CLOSE);
        // drawSurfaceNormals();
        pop()
    }

    drawSurfaceNormals(){
        // Draw surface normals for each line
        push()
        stroke('red');
        strokeWeight(2);
        for (let line_ of this.lines) {
            line_.drawNormal();
        }
        pop();
    }

    collides_with(mass_point_position){
        // check if the mass point is inside the bounding box
        if (mass_point_position.x > this.max_x || mass_point_position.x < this.min_x || mass_point_position.y > this.max_y || mass_point_position.y < this.min_y)
            return false;
        // check if the mass point is inside the polygon
        return this.surrounds(mass_point_position);
    }

    // checks is the position lies completely inside the polygon
    surrounds(position) {
        let count = 0;
        for (let line of this.lines){
            if (line.min_y > position.y || position.y > line.max_y)
                continue;
            
            // horizontal line that hits the mass point also intersects the polygon line
            // calculate the x of the intersection point
            // x - x1 = (x2 - x1) * (y - y1) / (y2 - y1) 
            let x = (line.q.x - line.p.x) * (position.y - line.p.y) / (line.q.y - line.p.y) + line.p.x;
            if (x < position.x)
                count++;
        }
        return count % 2 == 1;
    }

    collisionAlongRay(position, predictedPosition) {
        if (
            Math.max(position.x, predictedPosition.x) < this.min_x ||
            Math.min(position.x, predictedPosition.x) > this.max_x ||
            Math.max(position.y, predictedPosition.y) < this.min_y ||
            Math.min(position.y, predictedPosition.y) > this.max_y
        )
            return {
                collision: false,
                collidingLine: null,
                contact: null,
                t: Infinity,
            };
        
        // check if the ray from position to predictedPosition lies completely inside the polygon
        let pos_inside = this.surrounds(position);
        let pred_inside = this.surrounds(predictedPosition);
        if (pos_inside && pred_inside) {
            // continuous collision failed, do static collision
            let closest_line = this.get_closest_line(predictedPosition);
            let contact = closest_line.closest_point_on_line_from(predictedPosition);
            return {
                collision : true,
                collidingLine : closest_line,
                contact : contact,
                t: -1, // -1 indicates static collision
            }
        }
        
        let t = Infinity;
        let collision = false;
        let collidingLine = null;
        let contact = null;
        
        for (let line of this.lines) {
            let result = getLineIntersection(position, predictedPosition, line.p, line.q);
            if (!result.collision) continue;
            
            if (result.t < t){
                collision = result.collision;
                collidingLine = line;
                contact = result.contactPoint;
                t = result.t;
            }
        }
        if (t === Infinity) {
            return {
                collision: false,
                collidingLine: null,
                contact: null,
                t: Infinity,
            };
        } else {
            return {
                collision: true,
                collidingLine,
                contact,
                t: t,
            };
        }
    }

    get_closest_line(point){
        //! okay for line, not for line segment
        let min_d = Infinity;
        let closest_line;
        for (let line of this.lines) {
            let distance = line.distance_from_point(point)
            if (distance < min_d){
                min_d = distance;
                closest_line = line;
            }
        }
        return closest_line;
    }

    
    get_reflection_point_and_normal(predictedPosition, line){
        let closestPoint = line.closest_point_on_line_from(predictedPosition)
        let distance = p5.Vector.sub(closestPoint, predictedPosition).mag()
        let n_hat = closestLine.normal();
        return {
            point : closestPoint,
            dist: distance,
            normal : n_hat
        };
    }

    debugAlert(){
        fill('red');
        beginShape();
        for (let point of this.points){
            vertex(point.x, point.y);
        }
        endShape(CLOSE);
    }
}