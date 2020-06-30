#pragma once
#include <memory>
#include <vector>
#include <math.h>
#include <algorithm>

// -------------------------- //
// --- 2D ------------------- //
// -------------------------- //

namespace VPM
{

    struct IPoint2d
    {
        int x;
        int y;

        IPoint2d( ) {}

        IPoint2d( const int xy_ )
            :
                x( xy_ ),
                y( xy_ )
        {}

        IPoint2d( const int x_, const int y_ )
            :
                x( x_ ),
                y( y_ )
        {}
    };

    struct Point2d
    {
        double x;
        double y;

        Point2d( ) {}

        Point2d( const double xy_ )
            :
                x( xy_ ),
                y( xy_ )
        {}

        Point2d( const double x_, const double y_ )
            :
                x( x_ ),
                y( y_ )
        {}

        inline Point2d & operator+=(const Point2d& v )
        {
            x += v.x;
            y += v.y;
            return *this;
        }

        inline Point2d & operator-=(const Point2d& v )
        {
            x -= v.x;
            y -= v.y;
            return *this;
        }

    };

    //--------------//
    // Plus
    //--------------//
    inline Point2d operator+(const Point2d& lhs, const Point2d& rhs )
    {
        return Point2d(
                lhs.x + rhs.x,
                lhs.y + rhs.y
                );
    }

    inline Point2d operator+(const double a, const Point2d& v )
    {
        return Point2d(
                a + v.x,
                a + v.y
                );
    }

    inline Point2d operator+(const Point2d& v, const double a )
    {
        return a+v;
    }

    //--------------//
    // Minus
    //--------------//
    inline Point2d operator-(const double a, const Point2d& v )
    {
        return Point2d(
                a - v.x,
                a - v.y
                );
    }

    inline Point2d operator-(const Point2d& v, const double a )
    {
        return Point2d(
                v.x - a,
                v.y - a
                );
    }

    inline Point2d operator-(const Point2d& lhs, const Point2d& rhs )
    {
        return Point2d(
                lhs.x - rhs.x,
                lhs.y - rhs.y
                );
    }

    inline Point2d operator-(const Point2d& v)
    {
        return Point2d(-v.x, -v.y);
    }

    //--------------//
    // Times
    //--------------//
    inline Point2d operator*(const double a, const Point2d& rhs )
    {
        return Point2d(
                a*rhs.x,
                a*rhs.y
                );
    }
    inline Point2d operator*(const Point2d& lhs, const Point2d& rhs )
    {
        return Point2d(
                lhs.x*rhs.x,
                lhs.y*rhs.y
                );
    }

    inline Point2d operator*(const Point2d& v,const double a )
    {
        return a*v;
    }

    //--------------//
    // Devided
    //--------------//
    inline Point2d operator/(const double a, const Point2d& v )
    {
        return Point2d(
                a/v.x,
                a/v.y
                );
    }

    inline Point2d operator/(const Point2d& v,const double a )
    {
        return Point2d(
                v.x/a,
                v.y/a
                );
    }

    inline Point2d operator/(const Point2d& lhs, const Point2d& rhs )
    {
        return Point2d(
                lhs.x/rhs.x,
                lhs.y/rhs.y
                );
    }


    //--------------//
    // Algebraic operations
    //--------------//
    inline double dot(const Point2d& lhs, const Point2d& rhs )
    {
        return lhs.x*rhs.x + lhs.y*rhs.y;
    }

    inline Point2d perp(const Point2d & A)
    {
        return Point2d(A.y, -A.x);
    }

    //inline double cross(const Point2d A, const Point2d B)
    //{
    //    return A.x*B.y - A.y*B.x;
    //}

    inline double L2Norm(const Point2d A)
    {
        return sqrt(dot(A,A));
    }

}
