#pragma once
#include "Structure.hpp"

namespace VPM
{
    class Structure_Bar : public Structure
    {
    public:
       Structure_Bar();
        ~Structure_Bar();
        bool isInside(const Point3d pos, const double pad);

    private:

    };

}
