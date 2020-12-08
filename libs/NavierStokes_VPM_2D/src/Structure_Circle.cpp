#include "Structure_Circle.hpp"

namespace VPM
{

    Structure_Circle::Structure_Circle()
    {
       m_origo = Point2d(-.5, .111);
       m_radius = .2;
    }
    Structure_Circle::~Structure_Circle()
    {
    }

    bool Structure_Circle::isInside(const Point2d pos, const double pad)
    {

       bool ret = false;
       if (L2Norm(Point2d(pos-m_origo)) <= m_radius+pad)
       {
           ret = true;
       }
       return ret;
    }

    void Structure_Circle::getOrigo(Point2d & pos)
    {
        pos = m_origo;
    }

    double Structure_Circle::getCharacteristicLength()
    {
        return 2*m_radius;
    }


}
