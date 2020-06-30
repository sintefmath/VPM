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

    private:

    };

}
