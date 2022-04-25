#include "Structure_Transformed.hpp"

namespace VPM {
    Structure_Transformed::~Structure_Transformed() {}

    Structure_Transformed::Structure_Transformed(const Matrix2x2& linearTransformation, std::shared_ptr<Structure> structure) 
        : m_linearTransform(linearTransformation), m_structure(structure)
        {

        }
    
    bool Structure_Transformed::isInside(const Point2d pos, const double pad) {
        // First we transform the point
        auto transformedPoint = m_linearTransform * pos;
        return m_structure->isInside(transformedPoint, pad);

    }
}