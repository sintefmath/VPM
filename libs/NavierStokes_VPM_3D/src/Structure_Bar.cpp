#include "Structure_Bar.hpp"

namespace VPM
{

    Structure_Bar::Structure_Bar()
    {
    }
    Structure_Bar::~Structure_Bar()
    {
    }

    bool Structure_Bar::isInside(const Point3d pos, const double pad)
    {
       double circle_r = .2;
       Point3d circle_m = Point3d(-.5, .111, .011);

       bool ret = false;
       if (sqrt((pos.x-circle_m.x)*(pos.x-circle_m.x) + (pos.y-circle_m.y)*(pos.y-circle_m.y)) <= circle_r+pad)
       {
           ret = true;
       }
       return ret;
    }

}
