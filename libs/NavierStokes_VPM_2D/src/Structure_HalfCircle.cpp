#include "Structure_HalfCircle.hpp"

namespace VPM
{

    Structure_HalfCircle::Structure_HalfCircle()
    {
    }
    Structure_HalfCircle::~Structure_HalfCircle()
    {
    }

    bool Structure_HalfCircle::isInside(const Point2d pos, const double pad)
    {
       double circle_r = .2;
       Point2d circle_m = Point2d(-.5, .111);

       bool ret = false;
       if (L2Norm(Point2d(pos-circle_m)) <= circle_r+pad
           && pos.x<=circle_m.x+pad
          )
       {
           ret = true;
       }
       return ret;
    }

    void Structure_HalfCircle::getOrigo(Point2d & pos)
    {
    }

    double Structure_HalfCircle::getCharacteristicLength()
    {
        return 0.;// not correct
    }

}
