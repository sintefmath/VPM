#pragma once

#include "Particles3d.hpp"
#include "Parameters.hpp"

#include <vector>
#include <memory>

namespace ModulesHotel {
    class DiffusionEquationSolver3D;
}

namespace VPM
{
    class Split_Diffusion
    {
    public:
       Split_Diffusion();
        ~Split_Diffusion();


        void Euler_step(
                ParticleField & pf,
                const double dt
                );

        void RK2_step(
                ParticleField & pf,
                const double dt
                );

    private:

        int m_argc;
        char** m_argv;

        std::shared_ptr<ModulesHotel::DiffusionEquationSolver3D> m_diffeqsolver3d;

    };

}
