#include "Structure_Rectangle.hpp"

namespace VPM
{

    
    Structure_Rectangle::Structure_Rectangle(Point2d lower_corner, Point2d upper_corner)
    :
        m_lower_corner(lower_corner),
        m_upper_corner(upper_corner)
    {
    }
    Structure_Rectangle::~Structure_Rectangle()
    {
        // empty
    }

    bool Structure_Rectangle::isInside(const Point2d pos, const double pad)
    {

        return pos.x >= m_lower_corner.x && \
            pos.x <= m_upper_corner.x && \
            pos.y >= m_lower_corner.y && \
            pos.y <= m_upper_corner.y;
    }

}
