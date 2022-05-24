#include "Structure_Transformed.hpp"

namespace VPM {
    Structure_Transformed::~Structure_Transformed() {}

    Structure_Transformed::Structure_Transformed(const Matrix2x2& linearTransformation, const Point2d& translation, std::shared_ptr<Structure> structure) 
        : m_linearTransform(linearTransformation), m_translation(translation), m_structure(structure)
        {

        }
    
    bool Structure_Transformed::isInside(const Point2d pos, const double pad) {
        // First we transform the point
        auto transformedPoint = m_linearTransform * pos + m_translation;
        return m_structure->isInside(transformedPoint, pad);

    }
}