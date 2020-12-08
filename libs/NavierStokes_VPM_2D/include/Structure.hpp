#pragma once
#include "Point2d.hpp"

namespace VPM
{
    class Structure
    {
    public:
       Structure() {};
        ~Structure() {};
        virtual bool isInside(const Point2d pos, const double pad) = 0;
        virtual void getOrigo(Point2d & pos) = 0;
        virtual double getCharacteristicLength() = 0;

    private:

    };

}
