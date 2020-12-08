#include "Structure_Hilly.hpp"

namespace VPM
{

    Structure_Hilly::Structure_Hilly()
    {
    }
    Structure_Hilly::~Structure_Hilly()
    {
    }

    bool Structure_Hilly::isInside(const Point2d pos, const double pad)
    {
       bool ret = false;
       double xx = (pos.x-.25)/0.05;
       double yy = -.75 + .1*exp(-.1*xx*xx)+pad;
       if (pos.y<=yy)
       {
           ret = true;
       }
       return ret;
    }

    void Structure_Hilly::getOrigo(Point2d & pos)
    {
    }

    double Structure_Hilly::getCharacteristicLength()
    {
        return 0.;// not correct
    }

}
