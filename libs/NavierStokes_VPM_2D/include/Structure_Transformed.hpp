#pragma once
#include "Structure.hpp"
#include "Point2d.hpp"
#include "Matrix2x2.hpp"

namespace VPM
{
    class Structure_Transformed : public Structure
    {
    public:
        Structure_Transformed(const Matrix2x2& linearTransform, std::shared_ptr<Structure> structure);
        virtual ~Structure_Transformed();
        bool isInside(const Point2d pos, const double pad);

    private:
        Matrix2x2 m_linearTransform;        
        Point2d m_origo;
        std::shared_ptr<Structure> m_structure;
    };

}
