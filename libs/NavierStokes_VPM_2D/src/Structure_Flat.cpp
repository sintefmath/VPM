#include "Structure_Flat.hpp"

namespace VPM
{

    Structure_Flat::Structure_Flat()
    {
    }
    Structure_Flat::~Structure_Flat()
    {
    }

    bool Structure_Flat::isInside(const Point2d pos, const double pad)
    {
        bool ret = false;
        double yy = -.75+pad;
        if (pos.y<=yy)
        {
            ret = true;
        }
        return ret;

    }

    void Structure_Flat::getOrigo(Point2d & pos)
    {
    }

    double Structure_Flat::getCharacteristicLength()
    {
        return 0.;// not correct
    }

}
