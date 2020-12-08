#include "Structure_Ellipse.hpp"

namespace VPM
{

    Structure_Ellipse::Structure_Ellipse() :
        m_origo(Point2d(-.5,0)),
        m_semi_major_axis(.25),
        m_semi_major_axis2(m_semi_major_axis*m_semi_major_axis),
        m_semi_minor_axis(.15),
        m_semi_minor_axis2(m_semi_minor_axis*m_semi_minor_axis)
    {
    }
    Structure_Ellipse::Structure_Ellipse(Point2d origo, double semi_major_axis, double semi_minor_axis)
    :
        m_origo(origo),
        m_semi_major_axis(semi_major_axis),
        m_semi_major_axis2(semi_major_axis*semi_major_axis),
        m_semi_minor_axis(semi_minor_axis),
        m_semi_minor_axis2(semi_minor_axis*semi_minor_axis)
    {
    }
    Structure_Ellipse::~Structure_Ellipse()
    {
    }

    void Structure_Ellipse::getOrigo(Point2d & pos)
    {
        pos = m_origo;
    }

    double Structure_Ellipse::getCharacteristicLength()
    {
        return 0.2;// not correct
    }

    bool Structure_Ellipse::isInside(const Point2d pos, const double pad)
    {

        Point2d p = pos - m_origo;

        bool ret = false;

        double tmp = (p.x*p.x)/m_semi_major_axis2 + (p.y*p.y)/m_semi_minor_axis2;
        if (tmp <= 1 + pad)
        {
            ret = true;
        }

       return ret;
    }

}
