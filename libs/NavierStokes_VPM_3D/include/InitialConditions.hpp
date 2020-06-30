#pragma once
#include "Parameters.hpp"
#include "Particles3d.hpp"
#include <vector>

namespace VPM {

    //Point3d cartesian2polar(const Point3d x);
    //Point3d polar2cartesian(const Point3d x);

    Point3d omega_VortexRing_z(const Point3d x, const double strength=1., const double R0 = 1.5, const double eps0 = 0.3);
    Point3d omega_DoubleVortexRing_z(const Point3d x, const double strength=1.);

    Point3d omega_LambOseen_x(const Point3d x, const double strength=1., const double core_radius=0.15);
    Point3d omega_LambOseen_y(const Point3d x, const double strength=1., const double core_radius=0.15);
    Point3d omega_LambOseen_z(const Point3d x, const double strength=1., const double core_radius=0.15);
    Point3d omega_CounterRotatingVortexPair_x(const Point3d x, const double dist, const double strength=1., const double core_radius=0.15);
    Point3d omega_CounterRotatingVortexPair_y(const Point3d x, const double dist, const double strength=1., const double core_radius=0.15);
    Point3d omega_CounterRotatingVortexPair_z(const Point3d x, const double dist, const double strength=1., const double core_radius=0.15);

    Point3d omega_CoRotatingVortexPair_x(const Point3d x, const double dist, const double strength=1., const double core_radius=0.15);
    Point3d omega_CoRotatingVortexPair_y(const Point3d x, const double dist, const double strength=1., const double core_radius=0.15);
    Point3d omega_CoRotatingVortexPair_z(const Point3d x, const double dist, const double strength=1., const double core_radius=0.15);

    void init(
            const std::shared_ptr<VPM::Parameters> params,
            const int example_num,
            const double dist,
            const double strength,
            const double core_radius,
            std::vector<Point3d> & positions,
            std::vector<Point3d> & omega
            );

    void writeParticlesToFile(
            const std::string & filename,
            const int timestep,
            const ParticleField & pf
            );

    bool readParticlesFromFile(
            const std::string & filename,
            const std::shared_ptr<VPM::Parameters> params,
            std::vector<Point3d> & positions,
            std::vector<Point3d> & omega
            );
}

