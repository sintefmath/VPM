#pragma once
#include "Particles3d.hpp"

namespace VPM
{
    class Structure
    {
    public:
       Structure() {};
        ~Structure() {};
        virtual bool isInside(const Point3d pos, const double pad) = 0;

    private:

    };

}
