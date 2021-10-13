#include "Split_Diffusion.hpp"

#include "DiffusionEquationSolver2D.hpp"

namespace VPM
{

    Split_Diffusion::Split_Diffusion()
    {
        DiffusionEquationSolver2D_BC bc = DiffusionEquationSolver2D_BC(
                DIFFEQ2D_BC_NEUMANN,
                DIFFEQ2D_BC_NEUMANN,
                DIFFEQ2D_BC_NEUMANN,
                DIFFEQ2D_BC_NEUMANN
                );
        m_diffeqsolver2d = std::make_shared<ModulesHotel::DiffusionEquationSolver2D>(bc);

    };

    Split_Diffusion::~Split_Diffusion()
    {
        m_diffeqsolver2d.reset();
    };

    void Split_Diffusion::Euler_step(
            ParticleField & pf,
            const double dt
            )
    {

        if (pf.params.m_nu<1e-32)
        {
            return;
        }

        if (!pf.cartesianGrid)
        {
            std::cerr<<"  ---> Split_Diffusion:: Euler_step called on IRregular grid\n";
            exit(0);
            return;
        }

        #ifdef VPM_VERBOSE
        std::cerr<<"  ---> Split_Diffusion:: Euler_step\n";
        #endif

        m_diffeqsolver2d->solve(pf.omega,pf.params.m_num_px,pf.params.m_num_py, pf.params.m_dx, pf.params.m_dy, pf.params.m_nu, dt);

        pf.velocity_correspondsTo_omega = false;

    }

    void Split_Diffusion::RK2_step(
            ParticleField & pf,
            const double dt
            )
    {
        return;
    }

}
