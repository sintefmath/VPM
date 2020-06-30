#pragma once
#include <memory>
#include <vector>
#include <math.h>
#include <algorithm>
#include "Point3d.hpp"

// -------------------------- //
// --- 3D ------------------- //
// -------------------------- //

namespace VPM
{

    struct IPoint3d
    {
        int x;
        int y;
        int z;

        IPoint3d( ) {}

        IPoint3d( const int xyz_ )
            :
                x( xyz_ ),
                y( xyz_ ),
                z( xyz_ )
        {}

        IPoint3d( const int x_, const int y_, const int z_ )
            :
                x( x_ ),
                y( y_ ),
                z( z_ )
        {}
    };

    struct Point3d
    {
        double x;
        double y;
        double z;

        Point3d( ) {}

        Point3d( const double xyz_ )
            :
                x( xyz_ ),
                y( xyz_ ),
                z( xyz_ )
        {}

        Point3d( const double x_, const double y_, const double z_ )
            :
                x( x_ ),
                y( y_ ),
                z( z_ )
        {}

        inline Point3d & operator+=(const Point3d& v )
        {
            x += v.x;
            y += v.y;
            z += v.z;
            return *this;
        }

        inline Point3d & operator-=(const Point3d& v )
        {
            x -= v.x;
            y -= v.y;
            z -= v.z;
            return *this;
        }

    };

    //--------------//
    // Plus
    //--------------//
    inline Point3d operator+(const Point3d& lhs, const Point3d& rhs )
    {
        return Point3d(
                lhs.x + rhs.x,
                lhs.y + rhs.y,
                lhs.z + rhs.z
                );
    }

    inline Point3d operator+(const double a, const Point3d& v )
    {
        return Point3d(
                a + v.x,
                a + v.y,
                a + v.z
                );
    }

    inline Point3d operator+(const Point3d& v, const double a )
    {
        return a+v;
    }

    //--------------//
    // Minus
    //--------------//
    inline Point3d operator-(const double a, const Point3d& v )
    {
        return Point3d(
                a - v.x,
                a - v.y,
                a - v.z
                );
    }

    inline Point3d operator-(const Point3d& v, const double a )
    {
        return Point3d(
                v.x - a,
                v.y - a,
                v.z - a
                );
    }

    inline Point3d operator-(const Point3d& lhs, const Point3d& rhs )
    {
        return Point3d(
                lhs.x - rhs.x,
                lhs.y - rhs.y,
                lhs.z - rhs.z
                );
    }

    inline Point3d operator-(const Point3d& v)
    {
        return Point3d(-v.x, -v.y, -v.z);
    }

    //--------------//
    // Times
    //--------------//
    inline Point3d operator*(const double a, const Point3d& rhs )
    {
        return Point3d(
                a*rhs.x,
                a*rhs.y,
                a*rhs.z
                );
    }
    inline Point3d operator*(const Point3d& lhs, const Point3d& rhs )
    {
        return Point3d(
                lhs.x*rhs.x,
                lhs.y*rhs.y,
                lhs.z*rhs.z
                );
    }

    inline Point3d operator*(const Point3d& v,const double a )
    {
        return a*v;
    }

    //--------------//
    // Devided
    //--------------//
    inline Point3d operator/(const double a, const Point3d& v )
    {
        return Point3d(
                a/v.x,
                a/v.y,
                a/v.z
                );
    }

    inline Point3d operator/(const Point3d& v,const double a )
    {
        return Point3d(
                v.x/a,
                v.y/a,
                v.z/a
                );
    }

    inline Point3d operator/(const Point3d& lhs, const Point3d& rhs )
    {
        return Point3d(
                lhs.x/rhs.x,
                lhs.y/rhs.y,
                lhs.z/rhs.z
                );
    }


    //--------------//
    // Algebraic operations
    //--------------//
    inline double dot(const Point3d& lhs, const Point3d& rhs )
    {
        return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z;
    }

    //inline Point3d perp(const Point3d & A)
    //{
    //    return Point3d(A.y, -A.x);
    //}

    //inline double cross(const Point3d A, const Point3d B)
    //{
    //    return A.x*B.y - A.y*B.x;
    //}

    inline double L2Norm(const Point3d A)
    {
        return sqrt(dot(A,A));
    }

}
