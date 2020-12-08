#pragma once
#include "Structure.hpp"

namespace VPM
{
    class Structure_HalfCircle : public Structure
    {
    public:
       Structure_HalfCircle();
        ~Structure_HalfCircle();
        bool isInside(const Point2d pos, const double pad);
        void getOrigo(Point2d & pos);
        double getCharacteristicLength();

    private:

    };

}
