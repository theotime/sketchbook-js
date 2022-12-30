/*
    Physics
*/

class Physics {

    constructor() { }

    /*
        Construct an Align-Axis Bounding-Box in 3D.

        points = [x0, y0, z0, x1, y1, z1, ...]
        return [[x, y, z], [x, y, z]]
    */
    static computeAABB(points) {

        const numPoints = points.length / 3;

        let min = [points[0], points[1], points[2]];
        let max = [points[0], points[1], points[2]];

        for (let i = 1; i < numPoints; i++) {
            min[0] = Math.min(min[0], points[3 * i + 0]);
            min[1] = Math.min(min[1], points[3 * i + 1]);
            min[2] = Math.min(min[2], points[3 * i + 2]);

            max[0] = Math.max(max[0], points[3 * i + 0]);
            max[1] = Math.max(max[1], points[3 * i + 1]);
            max[2] = Math.max(max[2], points[3 * i + 2]);
        }

        let aabb = [];
        aabb[0] = min;
        aabb[1] = max;

        return aabb;
    }

    /*
        Line AABB intersection in 3D.
        Ref: https://tavianator.com/2011/ray_box.html

        origin, direction: [x, y, z]
        aabb: [[x, y, z], [x, y, z]]
        return [tmin, tmax]
    */
    static intersectLineAABB(origin, direction, aabb, result) {
        let tmin = -Infinity, tmax = Infinity;

        for (let i = 0; i < 3; i++) {
            let t0 = (aabb[0][i] - origin[i]) / direction[i]; // if /0 => Infinity
            let t1 = (aabb[1][i] - origin[i]) / direction[i];

            tmin = Math.max(tmin, Math.min(t0, t1));
            tmax = Math.min(tmax, Math.max(t0, t1));
        }

        result[0] = tmin; result[1] = tmax;
        return tmax >= tmin ? true : false;
    }
}