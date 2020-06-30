#include "Structure_InverseCircle.hpp"

namespace VPM
{

    Structure_InverseCircle::Structure_InverseCircle()
    {
    }
    Structure_InverseCircle::~Structure_InverseCircle()
    {
    }

    bool Structure_InverseCircle::isInside(const Point2d pos, const double pad)
    {
       double circle_r = .2;
       Point2d circle_m = Point2d(-.5, .111);

       bool ret = false;

        double dist = L2Norm(Point2d(pos-circle_m));
        if (pad>1e-6)
        {
            if (dist >= circle_r-pad && dist <= circle_r+pad)
            {
                ret = true;
            }
        }
        else
        {
            if (dist>=circle_r)
            {
                ret = true;
            }
        }

       return ret;
    }

}
