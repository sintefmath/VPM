#pragma once
#include "Parameters.hpp"
#include "Particles2d.hpp"
#include <vector>

namespace VPM {

    double omega_LambOseen(const Point2d x, const double strength=1., const double core_radius=0.15);
    double omega_CoRotatingVortexPair(const Point2d x, const double dist, const double strength=1., const double core_radius=0.15);
    double omega_CounterRotatingVortexPair(const Point2d x, const double dist, const double strength=1., const double core_radius=0.15);

    void init(
            const VPM::Parameters& params,
            const int example_num,
            const double dist,
            const double strength,
            const double core_radius,
            std::vector<Point2d> & positions,
            std::vector<double> & omega
            );

    void writeParticlesToFile(
            const std::string & filename,
            const int timestep,
            const ParticleField & pf
            );

    bool readParticlesFromFile(
            const std::string & filename,
            ParticleField & pf,
            bool random_velocity_dist = false
            );
}

