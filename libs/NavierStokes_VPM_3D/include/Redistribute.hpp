#pragma once
#include "Parameters.hpp"
#include "Particles3d.hpp"

namespace VPM
{
    enum REDISTRIBUTE3D { RED3D_G2P_OFF, RED3D_G2P_BUILD, RED3D_G2P_USE };
    class Redistribute
    {
    public:
        Redistribute();
        ~Redistribute();

        /*
         * changes pf.omega
         */
        void redistribute(
                ParticleField & pf,
                const REDISTRIBUTE3D g2p = RED3D_G2P_OFF
                );
        /*
         */
        template <class T>
        void redistribute(
                const std::vector<Point3d> & old_pos,
                const std::vector<Point3d> & new_pos,
                std::shared_ptr<Parameters> & params,
                std::vector<T> & field,
                const REDISTRIBUTE3D g2p = RED3D_G2P_OFF
                );
        /*
         */
        void get_indices(
                const std::shared_ptr<Parameters> & params,
                const int stenc,
                const VPM::Point3d & pos,
                VPM::IPoint3d & ind_min,
                VPM::IPoint3d & ind_max
                );
    private:

        /* lists particle numbers, e.g.:
         *     0,100,111,...        // grid point 0,0.0
         *     33,10,211,...        // grid point 1,0,0
         *     ...
         *     100,...        // grid point Nx, Ny, Nz
         */
        std::vector<std::vector<int>> m_grid_to_particles;

    };
}

