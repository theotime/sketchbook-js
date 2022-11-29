/*
    Matrix 4x4 and Vec3 operations.
    Should be enough for a small 3D engine. 
*/

// helpers for linear algebra basics
const Vec3 = {

    back: [0, 0, -1],
    down: [0, -1, 0],
    forward: [0, 0, 1],
    left: [-1, 0, 0],
    one: [1, 1, 1],
    right: [1, 0, 0],
    up: [0, 1, 0],
    zero: [0, 0, 0],

    copy: (v) => [...v],
    neg: (v) => [-v[0], -v[1], -v[2]],
    add: (v1, v2) => [v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]],
    sub: (v1, v2) => [v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]],
    scale: (v, s) => [v[0] * s, v[1] * s, v[2] * s],
    dot: (v1, v2) => v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2],
    length: (v) => Math.sqrt(Vec3.dot(v, v)),
    normalize: (v) => {
        let len = Vec3.length(v);
        return [v[0] / len, v[1] / len, v[2] / len]
    },
    // https://en.wikipedia.org/wiki/Cross_product#Coordinate_notation
    // a2*b3 - a3*b2, a3*b1 - a1*b3, a1*b2 - a2*b1
    cross: (v1, v2) => [
        v1[1] * v2[2] - v1[2] * v2[1],
        v1[2] * v2[0] - v1[0] * v2[2],
        v1[0] * v2[1] - v1[1] * v2[0]
    ]
};

class Mat4x4 {

    // column-major: computation & array declaration (OpenGL format)
    // https://www.scratchapixel.com/lessons/mathematics-physics-for-computer-graphics/geometry/row-major-vs-column-major-vector
    constructor(params = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]) { // default = identity
        // simple checks to catch basic mistakes
        if (!Array.isArray(params)) {
            throw "Not an array";
        } else if (params.some(isNaN)) {
            throw "Not a number.";
        } else if (params.length != 16) {
            throw "Not the right amount of params.";
        }

        this.value = params;
    }

    get() {
        return [...this.value];
    }

    copy() {
        return new Mat4x4(this.get());
    }

    getTranspose() {
        return this.copy().transpose();
    }

    getInverse(flag = false) {
        return this.copy().inverse(flag);
    }

    // https://developer.mozilla.org/en-US/docs/Web/API/WebGL_API/Matrix_math_for_the_web
    // converted to column-major
    multiplyVec4(v) {
        let x = this.value[0] * v[0] + this.value[4] * v[1] + this.value[8] * v[2] + this.value[12] * v[3];
        let y = this.value[1] * v[0] + this.value[5] * v[1] + this.value[9] * v[2] + this.value[13] * v[3];
        let z = this.value[2] * v[0] + this.value[6] * v[1] + this.value[10] * v[2] + this.value[14] * v[3];
        let w = this.value[3] * v[0] + this.value[7] * v[1] + this.value[11] * v[2] + this.value[15] * v[3];
        return [x, y, z, w];
    }

    multiplyMat4x4(matrix) {

        let m = matrix.get();

        let column0 = [m[0], m[1], m[2], m[3]];
        let column1 = [m[4], m[5], m[6], m[7]];
        let column2 = [m[8], m[9], m[10], m[11]];
        let column3 = [m[12], m[13], m[14], m[15]];

        let result0 = this.multiplyVec4(column0);
        let result1 = this.multiplyVec4(column1);
        let result2 = this.multiplyVec4(column2);
        let result3 = this.multiplyVec4(column3);

        let result = new Mat4x4([
            result0[0], result0[1], result0[2], result0[3],
            result1[0], result1[1], result1[2], result1[3],
            result2[0], result2[1], result2[2], result2[3],
            result3[0], result3[1], result3[2], result3[3]
        ]);

        return result;
    }

    transpose() {
        const swap = (a, m, n) => [a[m], a[n]] = [a[n], a[m]];

        swap(this.value, 1, 4);
        swap(this.value, 2, 8);
        swap(this.value, 3, 12);
        swap(this.value, 6, 9);
        swap(this.value, 7, 13);
        swap(this.value, 11, 14);

        return this;
    }

    /*
        Inverse of a 3x3 matrix (because it's simpler and enough for what we need now).

        XXX TODO: Handle affine transformations for Mat4x4 before the general case.

        Cramer's rule - Adjugate matrix
        https://gamemath.com/book/matrixmore.html
        https://en.wikipedia.org/wiki/Adjugate_matrix

        A * adj(A) = det(A) * I

        mat4 to mat3
        m0 m3 m6      m0 m4 m8  m12
        m1 m4 m7  <-  m1 m5 m9  m13
        m2 m5 m8      m2 m6 m10 m14
                      m3 m7 m11 m15

        m0 m3 m6   +  -  +
        m1 m4 m7   -  +  -
        m2 m5 m8   +  -  +

        Adjugate matrix is the transpose of the cofactor matrix.
        Cofactors:
        c0 = (+) m4m8 - m7m5
        c6 = (+) m1m5 - m4m2
        c2 = (+) m3m7 - m6m4
        c8 = (+) m0m4 - m3m1
        c4 = (+) m0m8 - m6m2
        c3 = (-) m1m8 - m7m2
        c7 = (-) m0m5 - m3m2
        c5 = (-) m0m7 - m6m1
        c1 = (-) m3m8 - m6m5

        det(A) = Sum [m(i,j) * c(i,j)]
        = m0 * c0 + m3 * c3 + m6 * c6
    */
    static BYPASS_EXCEPTION = true;
    inverse(bypassException = false) {
        // will throw an exception if matrix is not a 3x3 matrix into a 4x4 matrix.
        let is3x3matrix = false;
        if (this.value[3] == 0 && this.value[7] == 0 && this.value[11] == 0 &&
            this.value[12] == 0 && this.value[13] == 0 && this.value[14] == 0 &&
            this.value[15] == 1) {
            is3x3matrix = true;
        }

        if (!is3x3matrix && bypassException) {
            console.warn("Notice: Mat4x4::inverse() not fully implemented, will reset the 4th column & row with (0,0,0,1)");
        }
        if (!is3x3matrix && !bypassException) {
            console.warn(this.toString(2, "console"));
            throw "Error: Mat4x4::inverse() not fully implemented, 4th column & row must be (0,0,0,1)"
        }

        // convert to mat3 and not directly assigning into the mat4,
        // so we can refactor this function into Mat3x3.inverse() later when we'll fully implement Mat4x4.inverse().
        let m = [
            this.value[0], this.value[1], this.value[2], // column 1
            this.value[4], this.value[5], this.value[6], // column 2
            this.value[8], this.value[9], this.value[10] // column 3
        ];

        // compute cofactors
        let c = [];
        c[0] = m[4] * m[8] - m[7] * m[5];
        c[1] = m[6] * m[5] - m[3] * m[8];
        c[2] = m[3] * m[7] - m[6] * m[4];
        c[3] = m[7] * m[2] - m[1] * m[8];
        c[4] = m[0] * m[8] - m[7] * m[2];
        c[5] = m[6] * m[1] - m[0] * m[7];
        c[6] = m[1] * m[5] - m[4] * m[2];
        c[7] = m[3] * m[2] - m[0] * m[5];
        c[8] = m[0] * m[4] - m[3] * m[1];

        // compute inverse determinant
        let invDet = 1.0 / (m[0] * c[0] + m[3] * c[3] + m[6] * c[6]);

        // convert back to mat4, transpose to get adjugate matrix & compute the inverse in one pass
        // A^-1 = adj(A) / det(A)
        this.value = [
            c[0] * invDet, c[1] * invDet, c[2] * invDet, 0,
            c[3] * invDet, c[4] * invDet, c[5] * invDet, 0,
            c[6] * invDet, c[7] * invDet, c[8] * invDet, 0,
            0, 0, 0, 1,
        ];
        this.transpose(); // same as doing before the division

        return this;
    }

    translate(tx, ty, tz) {
        let m = new Mat4x4();
        this.value = m.fromTranslate(tx, ty, tz).multiplyMat4x4(this).get();
        return this;
    }

    rotate(angle, axis) {
        let m = new Mat4x4();
        this.value = m.fromAngleAxis(angle, axis).multiplyMat4x4(this).get();
        return this;
    }

    scale(sx, sy, sz) {
        let m = new Mat4x4();
        this.value = m.fromScale(sx, sy, sz).multiplyMat4x4(this).get();
        return this;
    }

    fromIdentity() {
        this.value = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
        return this;
    }

    // http://www.songho.ca/opengl/gl_quaternion.html
    //
    // 1-2yy-2zz     2xy-2sz      2xz+2sy      0
    // 2xy+2sz       1-2xx-2zz    2yz-2sx      0
    // 2xz-2sy       2yz+2sx      1-2xx-2yy    0
    // 0             0            0            1
    fromQuaternion(x, y, z, w) {
        // vector [x, y, z], scalar w
        this.value = [
            1 - 2 * y * y - 2 * z * z, 2 * x * y + 2 * w * z, 2 * x * z - 2 * w * y, 0, // column 1
            2 * x * y - 2 * w * z, 1 - 2 * x * x - 2 * z * z, 2 * y * z + 2 * w * x, 0, // column 2
            2 * x * z + 2 * w * y, 2 * y * z - 2 * w * x, 1 - 2 * x * x - 2 * y * y, 0, // column 3
            0, 0, 0, 1
        ];

        return this;
    }

    fromAngleAxis(angle, axis) {
        let cos = Math.cos(angle / 2.0);
        let sin = Math.sin(angle / 2.0);
        this.fromQuaternion(axis[0] * sin, axis[1] * sin, axis[2] * sin, cos);

        return this;
    }

    // 1 0 0 tx
    // 0 1 0 ty
    // 0 0 1 tz
    // 0 0 0 1
    fromTranslate(tx, ty, tz) {
        this.value = [
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            tx, ty, tz, 1
        ]; // column-major

        return this;
    }

    // http://www.songho.ca/opengl/gl_camera.html
    //
    // lx   ly   lz   -lx*ex-ly*ey-lz*ez
    // ux   uy   uz   -ux*ex-uy*ey-uz*ez
    // fx   fy   fz   -fx*ex-fy*ey-fz*ez
    // 0    0    0    1
    lookAt(eye, center, up) { // vec3
        let f = Vec3.normalize(Vec3.sub(eye, center));
        let l = Vec3.normalize(Vec3.cross(up, f));
        let u = Vec3.cross(f, l);

        this.value = [
            l[0], u[0], f[0], 0,
            l[1], u[1], f[1], 0,
            l[2], u[2], f[2], 0,
            Vec3.dot(l, Vec3.neg(eye)), Vec3.dot(u, Vec3.neg(eye)), Vec3.dot(f, Vec3.neg(eye)), 1
        ]; // column-major

        return this;
    }

    // http://www.songho.ca/opengl/gl_projectionmatrix.html
    // (with t = tanfov * n, r = t * aspect)
    //
    // S1   0   0   0
    // 0   S2   0   0
    // 0     0   A   B
    // 0     0   -1  0
    // S1 = 1/(aspect * tanfov)   S2 = 1/tanfov
    // A = -(f+n)/(f-n)   B = -(2*f*n)/(f-n)
    perspective(fov, aspect, near, far) {
        let tanfov = Math.tan((fov * Math.PI) / 360.0);
        let S1 = 1.0 / (aspect * tanfov);
        let S2 = 1.0 / tanfov;
        let A = -(far + near) / (far - near);
        let B = -(2 * far * near) / (far - near);
        this.value = [
            S1, 0, 0, 0,
            0, S2, 0, 0,
            0, 0, A, -1,
            0, 0, B, 0
        ];

        return this;
    }

    // sx  0  0   0
    // 0  sy  0   0
    // 0  0   sz  0
    // 0  0   0   1
    fromScale(sx, sy, sz) {
        this.value = [sx, 0, 0, 0, 0, sy, 0, 0, 0, 0, sz, 0, 0, 0, 0, 1];

        return this;
    }

    // XXX will work on a better implementation later
    almostEqual(matrix, threshold = Number.EPSILON) {
        for (let i = 0; i < 16; i++) {
            if (!(Math.abs(this.value[i] - matrix.value[i]) <= threshold)) {
                return false;
            }
        }
        return true;
    }

    round(rounding) {
        for (let i = 0; i < this.value.length; i++) {
            this.value[i] = this.value[i].toFixed(rounding);
            if (this.value[i] == 0) {
                this.value[i] = 0;
            } else {
                this.value[i] = parseFloat(this.value[i]);
            }
        }

        return this;
    }

    toString(rounding = 2, formatted = "html") {
        let string = "";

        let m = this.copy();
        m.round(rounding);

        let endl = " ";
        if (formatted == "html") {
            endl = "<br>";
        } else if (formatted == "console") {
            endl = "\n";
        }

        // first row
        string += m.value[0].toFixed(rounding) + " ";
        string += m.value[4].toFixed(rounding) + " ";
        string += m.value[8].toFixed(rounding) + " ";
        string += m.value[12].toFixed(rounding) + endl;

        // second row
        string += m.value[1].toFixed(rounding) + " ";
        string += m.value[5].toFixed(rounding) + " ";
        string += m.value[9].toFixed(rounding) + " ";
        string += m.value[13].toFixed(rounding) + endl;

        // third row
        string += m.value[2].toFixed(rounding) + " ";
        string += m.value[6].toFixed(rounding) + " ";
        string += m.value[10].toFixed(rounding) + " ";
        string += m.value[14].toFixed(rounding) + endl;

        // fourth row
        string += m.value[3].toFixed(rounding) + " ";
        string += m.value[7].toFixed(rounding) + " ";
        string += m.value[11].toFixed(rounding) + " ";
        string += m.value[15].toFixed(rounding) + endl;

        return string;
    }
}