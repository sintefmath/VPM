#pragma once
#include <memory>
#include <vector>
#include <math.h>
#include <algorithm>
#include "Point3d.hpp"
#include "Parameters.hpp"

// -------------------------- //
// --- 3D ------------------- //
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
        std::vector<Point3d> positions;
        std::vector<Point3d> regular_positions;
        std::vector<Point3d> omega;
        std::shared_ptr<Parameters> params;
        std::vector<Point3d> velocity;
        std::vector<Point3d> Jvel_x_omega;
        bool cartesianGrid;
        bool velocity_correspondsTo_omega;
        double Linfty_gradVelocity;
        double Linfty_omega;
    };

}
