#include "Point2d.hpp"

namespace VPM {
// TODO: Make this into a generic NxM minimatrix
struct Matrix2x2 {

    // Creates
    //   | a   b |
    //   | c   d |
    //
    Matrix2x2(double a, double b, double c, double d) 
        : a(a), b(b), c(c), d(d) {}
    
    Point2d operator*(const Point2d& point) {
        auto x = point.x;
        auto y = point.y;
        return Point2d( a * x + b * y, c * x + d * y );
    }

    static Matrix2x2 makeRotation(double angle) {
        return Matrix2x2(cos(angle), -sin(angle), sin(angle), cos(angle));
    }

    double a, b, c, d;

};
}