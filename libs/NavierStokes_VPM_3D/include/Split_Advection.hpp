#pragma once

#include "Particles3d.hpp"
#include "Parameters.hpp"
#include "Redistribute.hpp"
#include "Structure.hpp"

#include <vector>
#include <memory>

namespace ModulesHotel {
    class PoissonSolver3D;
}

namespace VPM
{
    class Split_Advection
    {
    public:
        Split_Advection(double lambda = 1.e4);

        ~Split_Advection();

        void Euler_step(
                ParticleField & pf,
                const double dt
                );

        void RK2_step(
                ParticleField & pf,
                const double dt
                );

        void setStructure(std::shared_ptr<VPM::Structure> structure);

        void getVelocityAndJacobian(
                ParticleField & pf,
                const double dt
                );


    private:
        int m_argc;
        char** m_argv;

        bool m_structure_set;

        std::shared_ptr<ModulesHotel::PoissonSolver3D> m_poissonsolver3d;

        std::shared_ptr<Redistribute> m_redistribute;
        std::shared_ptr<VPM::Structure> m_structure;

        void calcLinftyGradVelocity(
                ParticleField & pf
                );

        double m_lambda;


    };

}
