#pragma once
#include "Structure.hpp"

namespace VPM
{
    class Structure_Wing : public Structure
    {
    public:
       Structure_Wing();
        ~Structure_Wing();
        bool isInside(const Point2d pos, const double pad);
        void getOrigo(Point2d & pos);
        double getCharacteristicLength();

    private:
        Point2d m_origo;
        double m_charlength;
        std::vector<Point2d> m_p;

    };

}
