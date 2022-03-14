#include "UnionStructure.hpp"

namespace VPM {
    UnionStructure::UnionStructure(std::vector<std::shared_ptr<Structure>> structures) 
    : structures(structures)
    {
        // empty
    }

    bool UnionStructure::isInside(const Point2d pos, const double pad) {
        for (const auto& structure : structures) {
            if (structure->isInside(pos, pad)) {
                return true;
            }
        }
        return false;
    }

}