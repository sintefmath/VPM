#pragma once
#include "Structure.hpp"
#include "Point2d.hpp"

namespace VPM
{
    class Structure_Rectangle : public Structure
    {
    public:
        Structure_Rectangle(Point2d lower_corner, Point2d upper_corner);
        ~Structure_Rectangle();
        bool isInside(const Point2d pos, const double pad);

    private:
        Point2d m_lower_corner;
        Point2d m_upper_corner;
    };

}
