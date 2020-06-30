#pragma once
#include "Parameters.hpp"
#include "Particles2d.hpp"

namespace VPM
{
    enum REDISTRIBUTE2D { RED2D_G2P_OFF, RED2D_G2P_BUILD, RED2D_G2P_USE };
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
                const REDISTRIBUTE2D g2p = RED2D_G2P_OFF
                );
        /*
         */
        template <class T>
        void redistribute(
                const std::vector<Point2d> & old_pos,
                const std::vector<Point2d> & new_pos,
                Parameters & params,
                std::vector<T> & field,
                const REDISTRIBUTE2D g2p = RED2D_G2P_OFF
                );
        /*
         */
        void get_indices(
                const Parameters & params,
                const int stenc,
                const VPM::Point2d & pos,
                VPM::IPoint2d & ind_min,
                VPM::IPoint2d & ind_max
                );
    private:

        /* lists particle numbers, e.g.:
         *     0,100,111,...        // grid point 0,0
         *     33,10,211,...        // grid point 1,0
         *     ...
         *     100,...        // grid point Nx, Ny
         */
        std::vector<std::vector<int>> m_grid_to_particles;

    };
}

