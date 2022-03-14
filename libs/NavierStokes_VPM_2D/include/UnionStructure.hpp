#pragma once
#include "Structure.hpp"
#include <memory>
#include <vector>

namespace VPM
{
    class UnionStructure : public Structure
    {
    public:
        UnionStructure(std::vector<std::shared_ptr<Structure>> structures);

        bool isInside(const Point2d pos, const double pad);

    private:
        std::vector<std::shared_ptr<Structure>> structures;
    };

}
