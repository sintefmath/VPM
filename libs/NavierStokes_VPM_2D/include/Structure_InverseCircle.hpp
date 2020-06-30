#pragma once
#include "Structure.hpp"

namespace VPM
{
    class Structure_InverseCircle : public Structure
    {
    public:
       Structure_InverseCircle();
        ~Structure_InverseCircle();
        bool isInside(const Point2d pos, const double pad);

    private:

    };

}
