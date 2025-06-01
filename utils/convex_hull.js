function CCW(p0, p1, p2) {
    // returns 1 if p0, p1, p2 are in counter-clockwise order
    // returns -1 if they are in clockwise order
    // returns 0 if they are collinear
    let crossProduct = (p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x);
    return crossProduct > 0 ? -1 : (crossProduct < 0 ? 1 : 0); 
    // the inverse behaviour is because of the inverted y axis
  }

function getMaxYPoint(points) {
    // Find the point with the minimum y-coordinate
    // If there are multiple points with the same y-coordinate, return the one with the minimum x-coordinate
    let maxYPoint = points[0];
    for (let point of points) {
        if (point.y > maxYPoint.y || (point.y === maxYPoint.y && point.x < maxYPoint.x)) {
            maxYPoint = point;
        }
    }
    return maxYPoint;
}

function sortByAngle(points, lowest) {
    const distanceSq = (a, b) => (a.x - b.x) ** 2 + (a.y - b.y) ** 2;

    return points.slice().sort((a, b) => {
        if (a === lowest) return -1;
        if (b === lowest) return 1;

        const ccwResult = CCW(lowest, a, b);
        if (ccwResult != 0) return -ccwResult; // if ccw, then don't swap

        if (ccwResult === 0) {
            // Collinear: sort based on slope direction
            let dx = a.x - lowest.x;
            let dy = a.y - lowest.y;
            const m = dy / dx;
            if (m <= 0) {
                // Zero or positive slope → sort left to right (smaller x)
                return a.x - b.x;
            } else {
                // Negative or infinite slope → sort top to bottom (smaller y)
                return a.y - b.y;
            }
        }
    });
}

function grahamScan(points) {
    let maxYPoint = getMaxYPoint(points);
    let sorted = sortByAngle(points, maxYPoint);
    let stack = [];

    stack.push(sorted[0], sorted[1]);

    for (let i = 2; i < sorted.length; i++) {
        while (
            stack.length >= 2 &&
            CCW(stack[stack.length - 2], stack[stack.length - 1], sorted[i]) !== 1
        ) {
            stack.pop();
        }
        stack.push(sorted[i]);
    }

    return stack;
}