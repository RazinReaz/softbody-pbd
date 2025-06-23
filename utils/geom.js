function distSquared(p1, p2) {
    const dx = p1.x - p2.x;
    const dy = p1.y - p2.y;

    return dx ** 2 + dy ** 2;
}

function vecSetMult(result, vector, scalar) {
    result.set(vector.x, vector.y).mult(scalar);
}


function vecSetAdd(result, vec1, vec2) {
    result.set(vec1.x, vec1.y).add(vec2);
}

function vecSetSub(result, vec1, vec2) {
    result.set(vec1.x, vec1.y).sub(vec2);
}

