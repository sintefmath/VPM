#pragma once
#include "Structure.hpp"

namespace VPM
{
    class Structure_Hilly : public Structure
    {
    public:
       Structure_Hilly();
        ~Structure_Hilly();
        bool isInside(const Point2d pos, const double pad);

    private:

    };

}
