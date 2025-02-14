function polygon_points(x, y, n, radius){
    let points = [];
    let a = 2 * PI / n;
    for (let i = 0; i < n; i++){
        let r =  random(1, 2);
        points.push(createVector(x + r * radius * cos(i * a), y + r * radius * sin(i * a)));
    }
    return points;
}

class Polygon{
    constructor(x, y, n, radius = 50){
        this.points = [];
        if( Array.isArray(x) ){
            this.points = x;
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

        // // bounding box
        // push()
        // strokeWeight(3)
        // stroke('red')
        // line(this.min_x, this.min_y, this.min_x, this.max_y)
        // line(this.max_x, this.min_y, this.max_x, this.max_y)
        // line(this.min_x, this.min_y, this.max_x, this.min_y)
        // line(this.min_x, this.max_y, this.max_x, this.max_y)
        // pop()
    }

    collides_with(mass_point){
        // check if the mass point is inside the bounding box
        if (mass_point.position.x > this.max_x || mass_point.position.x < this.min_x || mass_point.position.y > this.max_y || mass_point.position.y < this.min_y)
            return false;
        // check if the mass point is inside the polygon
        let count = 0;
        for (let line of this.lines){
            if (line.min_y > mass_point.position.y || mass_point.position.y > line.max_y)
                continue;
            
            // horizontal line that hits the mass point also intersects the polygon line
            // calculate the x of the intersection point
            // x - x1 = (x2 - x1) * (y - y1) / (y2 - y1) 
            let x = (line.q.x - line.p.x) * (mass_point.position.y - line.p.y) / (line.q.y - line.p.y) + line.p.x;
            if (x < mass_point.position.x)
                count++;
        }
        return count % 2 == 1;
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

    
    get_reflection_point_and_normal(mass_point){
        let closest_line = this.get_closest_line(mass_point.position);

        return {
            point : closest_line.closest_point_on_line_from(mass_point.position),
            normal : closest_line.normal()
        };
    }

    debug_alert(){
        fill('red');
        beginShape();
        for (let point of this.points){
            vertex(point.x, point.y);
        }
        endShape(CLOSE);
    }
}