#include "Structure_Circle.hpp"

namespace VPM
{

    Structure_Circle::Structure_Circle()
    {
    }
    Structure_Circle::~Structure_Circle()
    {
    }

    bool Structure_Circle::isInside(const Point2d pos, const double pad)
    {
       double circle_r = .2;
       Point2d circle_m = Point2d(-.5, .111);

       bool ret = false;
       if (L2Norm(Point2d(pos-circle_m)) <= circle_r+pad)
       {
           ret = true;
       }
       return ret;
    }

}
