function isPrime(n) {
    if (n==1) return false;
    else if (n==2) return true;
    else if (n % 2 == 0) return false;
    else {
        for (let i = 3; i < Math.sqrt(n) + 1; i += 2) {
            if (n % i == 0)
                return false;
        }
        return true;
    }
}

function nextPrime(n) {
    for (let i = 0; ; i++) {
        if (isPrime(n+i))
            return n+i;
    }
}

function spatialHash(x, y) {
    const prime1 = 73856093;
    const prime2 = 19349663;
    return Math.abs((x * prime1) ^ (y * prime2));
}

function mortonHash(x, y) {
    return interleaveBits(x) | (interleaveBits(y) << 1);
}

function interleaveBits(n) {
    n = (n | (n << 8)) & 0x00FF00FF;
    n = (n | (n << 4)) & 0x0F0F0F0F;
    n = (n | (n << 2)) & 0x33333333;
    n = (n | (n << 1)) & 0x55555555;
    return n;
}


// for handling collision
class Grid {
    constructor(cellsize, numElements) {
        this.cellsize = cellsize;
        // this.columns = Math.ceil(screenWidth / cellsize);
        // this.rows = Math.ceil(screenHeight / cellsize);
        this.tablesize = nextPrime(2 * numElements)
  
        // the sparse solution
        this.cellsArray = [];
        for (let i = 0; i < this.tablesize; i++) {
            this.cellsArray.push([]);
        }
        this.size = 0;
  
  
        // // the dense solution
        // this.cellStart = new Array(this.tablesize + 1).fill(0);
        // this.cellEntries = new Array(numElements).fill(null);
    }
    
    #getCellIndices(elementPosition) {
        // const {x:vx, y:vy} = element.getVelocity();
        const xi = Math.floor(elementPosition.x / this.cellsize);
        const yi = Math.floor(elementPosition.y / this.cellsize);
        // console.log(`inside #getCellIndicies: xi ${xi} and yi ${yi}`)
        return {xi, yi};
    }

    #hashFunction(x, y) {
        return spatialHash(x, y);
    }
    
    #getArrayIndex(x, y) {
        let idx = this.#hashFunction(x, y) % this.tablesize;
        if (idx < 0){
            idx += this.tablesize;
        }

        return idx;
    }

    update(element){
        const {xi:new_xi, yi:new_yi} = this.#getCellIndices(element.position)
        const arrayIndex = this.#getArrayIndex(new_xi, new_yi);
        
        if (arrayIndex === element.cellKey) return;

        if (element.cellKey !== - 1){
            this.remove(element);
        }
        this.add(element, arrayIndex);

    }
  
    add(element, arrayIndex = -1) {
        if (arrayIndex !== -1) {
            const {xi, yi} = this.#getCellIndices(element.position)
            const arrayIndex = this.#getArrayIndex(xi, yi);
        }
        this.cellsArray[arrayIndex].push(element);
        element.cellKey = arrayIndex;
        this.size ++;
    }
  
    remove(element) {
        const arrayIndex = element.cellKey;
        let cell = this.cellsArray[arrayIndex];
        // cell is an array
        const elementIndex = cell.indexOf(element);

        if (elementIndex !== -1) {
            cell.splice(elementIndex, 1);
            this.size--;
        } else {
            console.warn("Grid: Tried to remove element not found in its cell.");
        }
    
        // element.cellKey = -1; //! reset?
    }

    query(elementPosition){
        let q = [];
        const {xi, yi} = this.#getCellIndices(elementPosition); //! maybe I can only pass the posiition?
        // console.log(`inside query: element.position: ${element.position.x}, ${element.position.y}`)
        // console.log(`inside query: xi ${xi} and yi ${yi}`)
        for (let x = xi - 1; x <= xi + 1; x++){
            for (let y = yi - 1; y <= yi + 1; y++) {
                let idx = this.#getArrayIndex(x, y);
                let neighbours = this.cellsArray[idx];
                // console.log(`inside query: neighbours length: ${neighbours.length}`)
                if (neighbours.length == 0) continue;
                q.push(...neighbours);
            }
        }
        q = q.filter(e => e.position !== elementPosition);
        return q;
    }
}