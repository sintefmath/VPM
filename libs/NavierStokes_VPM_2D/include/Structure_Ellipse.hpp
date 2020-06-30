#pragma once
#include "Structure.hpp"
#include "Point2d.hpp"

namespace VPM
{
    class Structure_Ellipse : public Structure
    {
    public:
        Structure_Ellipse();
        Structure_Ellipse(Point2d origo, double semi_major_axis, double semi_minor_axis);
        ~Structure_Ellipse();
        bool isInside(const Point2d pos, const double pad);

    private:

        double m_semi_major_axis;
        double m_semi_minor_axis;
        double m_semi_major_axis2;
        double m_semi_minor_axis2;
        Point2d m_origo;

    };

}
