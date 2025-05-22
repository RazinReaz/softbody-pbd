function polygon_points(x, y, n, radius){
    let points = [];
    let a = 2 * PI / n;
    for (let i = 0; i < n; i++){
        let r =  random(1, 2);
        points.push(createVector(x + r * radius * cos(i * a), y - r * radius * sin(i * a)));
    }
    return points;
}

// Orientation helper
function getOrientation(a, b, c) {
// a, b, c are all 2D vectors
    const crossProd = (b.y - a.y) * (c.x - b.x) - (b.x - a.x) * (c.y - b.y);
    if (crossProd === 0) return 0;  // colinear
    return crossProd > 0 ? 1 : 2;   // 1 = clockwise, 2 = counterclockwise
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
    if (denom == 0) return {collision: false, contactPoint: null}; // colinear

    let denomPositive = denom > 0;

    let s02 = {
        x: p0.x - p2.x,
        y: p0.y - p2.y
    }

    let s_numer = s10.x * s02.y - s10.y * s02.x;
    if ((s_numer < 0) == denomPositive)
        return {collision: false, contactPoint: null}; // no collision

    let t_numer = s32.x * s02.y - s32.y * s02.x;
    if ((t_numer < 0) === denomPositive)
        return { collision: false, contactPoint: null }; // No collision

    if (((s_numer > denom) === denomPositive) || ((t_numer > denom) === denomPositive))
        return { collision: false, contactPoint: null }; // No collision

    // Collision detected
    let t = t_numer / denom;
    let intersection = {
        x: p0.x + (t * s10.x),
        y: p0.y + (t * s10.y)
    };
    // console.log("intersection: ", intersection)
    // console.log("p0.x + (t * s10.x): ", p0.x + (t * s10.x))
    // console.log("p0.y + (t * s10.y): ", p0.y + (t * s10.y))

    return { collision: true, contactPoint: intersection };
}



class Polygon{
    constructor(x, y, n, radius = 50){
        this.points = [];
        if (Array.isArray(x)) {
            for (let pt of x) {
                this.points.push(pt)
            }
        } else {
            this.points = polygon_points(x, y, n, radius);
        }
        this.lines = [];
        for (let i = 0; i < this.points.length; i++){
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

    show(fillColor = 150){
        // TODO fill using other methods
        push()
        fill(fillColor)
        beginShape();
        for (let point of this.points){
            vertex(point.x, point.y);
        }
        endShape(CLOSE);
        pop()
    }

    collides_with(mass_point_position){
        // check if the mass point is inside the bounding box
        if (mass_point_position.x > this.max_x || mass_point_position.x < this.min_x || mass_point_position.y > this.max_y || mass_point_position.y < this.min_y)
            return false;
        // check if the mass point is inside the polygon
        let count = 0;
        for (let line of this.lines){
            if (line.min_y > mass_point_position.y || mass_point_position.y > line.max_y)
                continue;
            
            // horizontal line that hits the mass point also intersects the polygon line
            // calculate the x of the intersection point
            // x - x1 = (x2 - x1) * (y - y1) / (y2 - y1) 
            let x = (line.q.x - line.p.x) * (mass_point_position.y - line.p.y) / (line.q.y - line.p.y) + line.p.x;
            if (x < mass_point_position.x)
                count++;
        }
        return count % 2 == 1;
    }

    continuousCollisionAlong(position, predictedPosition) {
        if (
            Math.max(position.x, predictedPosition.x) < this.min_x ||
            Math.min(position.x, predictedPosition.x) > this.max_x ||
            Math.max(position.y, predictedPosition.y) < this.min_y ||
            Math.min(position.y, predictedPosition.y) > this.max_y
        )
            return {
                collision: false,
                collidingLine: null,
                contact: null
            };

        for (let line of this.lines) {
            let result = getLineIntersection(position, predictedPosition, line.p, line.q);
            if (!result.collision) continue;

            return {
                collision: result.collision,
                collidingLine: line,
                contact: result.contactPoint
            }
        //     const o1 = getOrientation(position, predictedPosition, line.p);
        //     const o2 = getOrientation(position, predictedPosition, line.q);
        //     const o3 = getOrientation(line.p, line.q, position);
        //     const o4 = getOrientation(line.p, line.q, predictedPosition);

        //     if (o1 !== o2 && o3 !== o4 ) 
        //         return {
        //             collision: true,
        //             collidingLine: line
        //         };

        //     if (o1 === 0 && onSegment(position, predictedPosition, line.p)) 
        //         return {
        //             collision: true,
        //             collidingLine: line
        //         };

        //     if (o2 === 0 && onSegment(position, predictedPosition, line.q)) 
        //         return {
        //             collision: true,
        //             collidingLine: line
        //         };

        //     if (o3 === 0 && onSegment(line.p, line.q, position)) 
        //         return {
        //             collision: true,
        //             collidingLine: line
        //         };

        //     if (o4 === 0 && onSegment(line.p, line.q, predictedPosition)) 
        //         return {
        //             collision: true,
        //             collidingLine: line
        //         };

        }
        return {
            collision: false,
            collidingLine: null,
            contact: null
        };
    }

    get_closest_line(point){
        //! okay for line, not for line segment
        let min_d = 10000;
        let closest_line;
        for (let line of this.lines) {
            if (line.distance_from_point(point) < min_d){
                min_d = line.distance_from_point(point);
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