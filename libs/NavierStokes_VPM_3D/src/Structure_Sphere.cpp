#include "Structure_Sphere.hpp"

namespace VPM
{

    Structure_Sphere::Structure_Sphere()
    {
    }
    Structure_Sphere::~Structure_Sphere()
    {
    }

    bool Structure_Sphere::isInside(const Point3d pos, const double pad)
    {
       double circle_r = .2;
       Point3d circle_m = Point3d(-.5, .011, .011);

       bool ret = false;
       if (L2Norm(Point3d(pos-circle_m)) <= circle_r+pad)
       {
           ret = true;
       }
       return ret;
    }

}
