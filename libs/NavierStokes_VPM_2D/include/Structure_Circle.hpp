#pragma once
#include "Structure.hpp"

namespace VPM
{
    class Structure_Circle : public Structure
    {
    public:
       Structure_Circle();
        ~Structure_Circle();
        bool isInside(const Point2d pos, const double pad);
        void getOrigo(Point2d & pos);
        double getCharacteristicLength();

    private:
        Point2d m_origo;
        double m_radius;

    };

}
