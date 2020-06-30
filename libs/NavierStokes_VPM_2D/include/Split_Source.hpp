#pragma once
#include "Parameters.hpp"
#include "Point2d.hpp"

#include <vector>
#include <memory>

namespace VPM
{
    class Split_Source
    {
    public:
       Split_Source(
               );
        ~Split_Source();
        std::vector<double> step(
                const std::vector<double> & omega,
                const std::shared_ptr<Parameters> params
                );

        void setVelocityFields(
                const std::vector<Point2d> & velocity,
                const std::vector<Point2d> & background_velocity,
                const std::shared_ptr<Parameters> params
                );

    private:

        //std::vector<Point2d> m_velocity;
        //std::vector<Point2d> m_background_velocity;
        std::vector<Point2d> m_total_velocity;
        std::vector<double> m_curl_background_velocity;

    };

}
