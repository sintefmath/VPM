#pragma once
#include <memory>
#include <vector>
#include <math.h>
#include <algorithm>
#include "Point2d.hpp"
#include "Parameters.hpp"

// -------------------------- //
// --- 2D ------------------- //
// -------------------------- //

namespace VPM
{

    /*
     * Important note:
     * In order to save CPU/memory operations, beware:
     * if (cartesianGrid)
     * {
     *       valid: regular_positions
     *     INvalid: positions
     * }
     * else
     * {
     *       valid: positions
     *     INvalid: regular_positions
     * }
     *    
     */
    struct ParticleField {
        std::vector<Point2d> positions;
        std::vector<Point2d> regular_positions;
        std::vector<double> omega;
        Parameters params;
        std::vector<Point2d> velocity;
        bool cartesianGrid;
        bool velocity_correspondsTo_omega;
        double Linfty_gradVelocity;
        double time;
    };

}
