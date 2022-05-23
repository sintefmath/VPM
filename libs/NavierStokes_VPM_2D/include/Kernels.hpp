#pragma once
#include "Parameters.hpp"
#include "Point2d.hpp"
#define PSEorder 2
// ugly hack...
double global_eps;
double global_delta;
double global_sigma;
#include "kernel_Base.hpp"
#include "H2_2D_Node.hpp"
typedef Point PosType;


inline PosType operator/(const double a, const PosType& v )
{
    return PosType(
        a/v.x,
        a/v.y
#ifdef _3D_VPM_
        ,
        a/v.z
#endif
        );
}
inline PosType operator/(const PosType& v,const double a )
{
    return PosType(
        v.x/a,
        v.y/a
#ifdef _3D_VPM_
        ,
        v.z/a
#endif
        );
}

#ifdef _3D_VPM_
inline PosType operator+(const PosType& lhs, const PosType& rhs )
{
    return PosType(
        lhs.x + rhs.x,
        lhs.y + rhs.y,
        lhs.z + rhs.z
        );
}
#endif

inline PosType operator+(const double a, const PosType& v )
{
    return PosType(
        a + v.x,
        a + v.y
#ifdef _3D_VPM_
        ,
        a + v.z
#endif
        );
}
inline PosType operator+(const PosType& v, const double a )
{
    return a+v;
}
inline PosType operator-(const double a, const PosType& v )
{
    return PosType(
        a - v.x,
        a - v.y
#ifdef _3D_VPM_
        ,
        a - v.z
#endif
        );
}
inline PosType operator-(const PosType& v, const double a )
{
    return PosType(
        v.x - a,
        v.y - a
#ifdef _3D_VPM_
        ,
        v.z - a
#endif
        );
}
inline PosType operator-(const PosType& lhs, const PosType& rhs )
{
    return PosType(
        lhs.x - rhs.x,
        lhs.y - rhs.y
#ifdef _3D_VPM_
        ,
        lhs.z - rhs.z
#endif
        );
}
inline PosType operator-(const PosType& v)
{
    return PosType(-v.x, -v.y
#ifdef _3D_VPM_
                   , -v.z
#endif
                   );
}
inline PosType operator*(const double a, const PosType& rhs )
{
    return PosType(
        a*rhs.x,
        a*rhs.y
#ifdef _3D_VPM_
        ,
        a*rhs.z
#endif
        );
}
inline PosType operator*(const PosType& v,const double a )
{
    return a*v;
}

inline double dot(const PosType & lhs, const PosType& rhs )
{
    return lhs.x*rhs.x + lhs.y*rhs.y;
}
inline double moll(const double r2, const double eps)
{
    // with r = x^2 + y^2
    // x, y |-> (1-( 1 - 2*(r/eps)^2 + 1/2*(r/eps)^4)*exp(-(r/eps)^2) )
    double r_eps2 = r2/(eps*eps);
    return 1. - (1. - 2.*r_eps2 + double(0.5)*r_eps2*r_eps2) * exp(-r_eps2);
}


inline double eta(const PosType & p, const double sigma)
{
    double r2 = dot(p,p)/sigma;
#if PSEorder==2
    // 2nd order:
    return 4./M_PI*exp(-r2);
#endif
#if PSEorder==4
    // 4th order:
    return 4./M_PI*(3-r2)*exp(-r2);
#endif
#if PSEorder==6
    // 6th order:
    return 2./M_PI*(12-8*r2+r2*r2)*exp(-r2);
#endif
#if PSEorder==8
    // 8th order:
    double r4 = r2*r2;
    double r6 = r4*r2;
    double r8 = r6*r2;
    return 1./(6*M_PI)*(360-480*r2+180*r4-24*r6+r8)*exp(-r2);
#endif
    //    //return 1./(4*M_PI)*exp(-r2/4.);
    //    if (r2<=4)
    //    {
    //        double C = 0.835;
    //        return C/(1.+r2);
    //    }
    //    else
    //    {
    //        return double(0.);
    //    }
}



/*
     * From: Vortex Methods: Theory and Practice, G. Cottet, and P. Koumoutsakos, Table 2.1 p. 27.
     */

class Kernel_K2_order6_x: public kernel_Base
{
public:
    virtual double kernel_Func(Point p0, Point p1) override
    {
        Point p = p0-p1;
        double r2 = dot(p,p);
        if (r2<1e-5)
        {
            return 0.;
        }
        return -p.y/(2*M_PI*r2)*moll(r2,global_eps);
    }
};
class Kernel_K2_order6_y: public kernel_Base
{
public:
    virtual double kernel_Func(Point p0, Point p1) override
    {
        Point p = p0-p1;
        double r2 = dot(p,p);
        if (r2<1e-5)
        {
            return 0.;
        }
        return +p.x/(2*M_PI*r2)*moll(r2,global_eps);;
    }
};
class Kernel_diffusion_eta: public kernel_Base
{
public:
    virtual double kernel_Func(Point p0, Point p1) override
    {
        return eta(p0-p1, global_sigma);
    }
};
