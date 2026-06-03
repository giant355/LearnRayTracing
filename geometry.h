#pragma once
#include <cmath>
#include <cassert>
#include <iostream>

template<int n> struct vec {
    double data[n] = { 0 };
    double& operator[](const int i) { assert(i >= 0 && i < n); return data[i]; }
    double  operator[](const int i)const { assert(i >= 0 && i < n); return data[i]; }
};

template<> struct vec<2> {
    double x = 0, y = 0;
    double& operator[](int idx) { assert(idx >= 0 && idx < 2); return idx == 0 ? x : y; };
    double  operator[](int idx)const { assert(idx >= 0 && idx < 2); return idx == 0 ? x : y; };
};

template<> struct vec<3> {
    double x = 0, y = 0, z = 0;
    double& operator[](int idx) { assert(idx >= 0 && idx < 3); return idx == 0 ? x : (idx == 1 ? y : z); }
    double  operator[](int idx)const { assert(idx >= 0 && idx < 3); return idx == 0 ? x : (idx == 1 ? y : z); }
};

template<> struct vec<4> {
    double x = 0, y = 0, z = 0, w = 0;
    double& operator[](const int idx) { assert(idx >= 0 && idx < 4); return idx == 0 ? x : (idx == 1 ? y : (idx == 2 ? z : w)); }
    double  operator[](const int idx)const { assert(idx >= 0 && idx < 4); return idx == 0 ? x : (idx == 1 ? y : (idx == 2 ? z : w)); }
    vec<2> xy() const { return { x, y }; }
    vec<3> xyz() const { return { x, y, z }; }
};

typedef vec<2> vec2;
typedef vec<3> vec3;
typedef vec<4> vec4;

template<int n> double operator*(const vec<n>& lhs, const vec<n>& rhs) {
    double ret = 0;
    for (int i = 0; i < n; ++i) ret += lhs[i] * rhs[i];
    return ret;
}

template<int n> vec<n> operator+(const vec<n>& lhs, const vec<n>& rhs) {
    vec<n> ret = lhs;
    for (int i = 0; i < n; ++i) ret[i] += rhs[i];
    return ret;
}

template<int n> vec<n> operator-(const vec<n>& lhs, const vec<n>& rhs) {
    vec<n> ret = lhs;
    for (int i = 0; i < n; ++i) ret[i] -= rhs[i];
    return ret;
}

template<int n> vec<n> operator*(const vec<n>& lhs, const double& rhs) {
    vec<n> ret = lhs;
    for (int i = 0; i < n; ++i) ret[i] *= rhs;
    return ret;
}

template<int n> vec<n> operator*(const double& lhs, const vec<n>& rhs) {
    return rhs * lhs;
}

template<int n> vec<n> operator/(const vec<n>& lhs, const double& rhs) {
    vec<n> ret = lhs;
    for (int i = 0; i < n; ++i) ret[i] /= rhs;
    return ret;
}

template<int n> std::ostream& operator<<(std::ostream& lhs, const vec<n>& rhs) {
    for (int i = 0; i < n; ++i) lhs << rhs[i] << " ";
    return lhs;
}

template<int n> double norm(const vec<n>& v) { return std::sqrt(v * v); }
template<int n> vec<n> normalized(const vec<n>& v) { return v / norm(v); }

inline vec3 cross(const vec3& v1, const vec3& v2) {
    return { v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x };
}

template<int nrows, int ncols> struct mat;
template<int ncols> struct dt;

template<int nrows, int ncols> struct mat {
    vec<ncols> rows[nrows] = { {} };
    vec<ncols>& operator[](int idx) { assert(idx >= 0 && idx < nrows); return rows[idx]; }
    const vec<ncols>& operator[](int idx)const { assert(idx >= 0 && idx < nrows); return rows[idx]; }

    double det()const {
        static_assert(nrows == ncols, "Matrix must be square to calculate determinant.");
        return dt<ncols>::det(*this);
    }

    double cofactor(int row, int col) const {
        mat<nrows - 1, ncols - 1> subMatrix;
        for (int i = 0; i < nrows - 1; ++i) {
            for (int j = 0; j < ncols - 1; ++j) {
                subMatrix[i][j] = rows[i + (i >= row)][j + (j >= col)];
            }
        }
        return ((row + col) % 2 ? -1 : 1) * subMatrix.det();
    }

    mat<nrows, ncols> invert_transpose()const {
        static_assert(nrows == ncols, "Matrix must be square to invert.");
        mat<nrows, ncols> cofactorMatrix;
        for (int i = 0; i < nrows; ++i)
            for (int j = 0; j < ncols; ++j)
                cofactorMatrix[i][j] = cofactor(i, j);
        return cofactorMatrix / det();
    }

    mat<nrows, ncols> invert()const {
        return invert_transpose().transpose();
    }

    mat<ncols, nrows> transpose()const {
        mat<ncols, nrows> ret;
        for (int i = 0; i < nrows; ++i)
            for (int j = 0; j < ncols; ++j)
                ret[j][i] = rows[i][j];
        return ret;
    }
};

template<int nrows, int ncols> vec<ncols> operator*(const vec<nrows>& lhs, const mat<nrows, ncols>& rhs) {
    vec<ncols> ret;
    for (int j = 0; j < ncols; ++j) {
        ret[j] = 0;
        for (int i = 0; i < nrows; ++i) ret[j] += lhs[i] * rhs[i][j];
    }
    return ret;
}

template<int nrows, int ncols> vec<nrows> operator*(const mat<nrows, ncols>& lhs, const vec<ncols>& rhs) {
    vec<nrows> ret;
    for (int i = 0; i < nrows; ++i) {
        ret[i] = 0;
        for (int j = 0; j < ncols; ++j) ret[i] += rhs[j] * lhs[i][j];
    }
    return ret;
}

template<int R1, int C1, int C2> mat<R1, C2> operator*(const mat<R1, C1>& lhs, const mat<C1, C2>& rhs) {
    mat<R1, C2> ret;
    for (int i = 0; i < R1; ++i) {
        for (int j = 0; j < C2; ++j) {
            ret[i][j] = 0;
            for (int k = 0; k < C1; ++k) ret[i][j] += lhs[i][k] * rhs[k][j];
        }
    }
    return ret;
}

template<int nrows, int ncols> mat<nrows, ncols> operator*(const mat<nrows, ncols>& lhs, const double& val) {
    mat<nrows, ncols> result;
    for (int i = 0; i < nrows; ++i) result[i] = lhs[i] * val;
    return result;
}

template<int nrows, int ncols> mat<nrows, ncols> operator/(const mat<nrows, ncols>& lhs, const double& val) {
    mat<nrows, ncols> result;
    for (int i = 0; i < nrows; ++i) result[i] = lhs[i] / val;
    return result;
}

template<int nrows, int ncols> std::ostream& operator<<(std::ostream& lhs, const mat<nrows, ncols>& rhs) {
    for (int i = 0; i < nrows; i++) lhs << rhs[i] << std::endl;
    return lhs;
}

template<int ncols> struct dt {
    static double det(const mat<ncols, ncols>& m) {
        double ret = 0;
        for (int i = 0; i < ncols; ++i) ret += m[0][i] * m.cofactor(0, i);
        return ret;
    }
};

template<> struct dt<1> {
    static double det(const mat<1, 1>& m) {
        return m[0][0];
    }
};