#pragma once
#include "Structure.hpp"

namespace VPM
{
    class Structure_Flat : public Structure
    {
    public:
       Structure_Flat();
        ~Structure_Flat();
        bool isInside(const Point3d pos, const double pad);

    private:

    };

}
