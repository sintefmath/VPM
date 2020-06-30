#pragma once
#include "Structure.hpp"

namespace VPM
{
    class Structure_Sphere : public Structure
    {
    public:
       Structure_Sphere();
        ~Structure_Sphere();
        bool isInside(const Point3d pos, const double pad);

    private:

    };

}
