#include "Split_Diffusion.hpp"

#include "DiffusionEquationSolver3D.hpp"

namespace VPM
{

    Split_Diffusion::Split_Diffusion()
    {
        DiffusionEquationSolver3D_BC bc = DiffusionEquationSolver3D_BC(
                DIFFEQ3D_BC_NEUMANN,
                DIFFEQ3D_BC_NEUMANN,
                DIFFEQ3D_BC_NEUMANN,
                DIFFEQ3D_BC_NEUMANN,
                DIFFEQ3D_BC_NEUMANN,
                DIFFEQ3D_BC_NEUMANN
                );
        m_diffeqsolver3d = std::make_shared<ModulesHotel::DiffusionEquationSolver3D>(bc);
    };

    Split_Diffusion::~Split_Diffusion()
    {
        m_diffeqsolver3d.reset();
    };


    void Split_Diffusion::Euler_step(
            ParticleField & pf,
            const double dt
            )
    {

        if (pf.params->m_nu<1e-32)
        {
            return;
        }

        if (!pf.cartesianGrid)
        {
            std::cerr<<"  ---> Split_Diffusion:: Euler_step called on IRregular grid\n";
            exit(0);
            return;
        }

        std::cerr<<"  ---> Split_Diffusion:: Euler_step\n";


        std::vector<double> tmp;
        tmp.clear();
        tmp.resize(pf.positions.size(), 0.);

        // x-dir
        for(unsigned int i=0; i<pf.positions.size(); i++ )
        {
           tmp[i] = pf.omega[i].x;
        }

        m_diffeqsolver3d->solve(tmp,pf.params->m_num_px,pf.params->m_num_py,pf.params->m_num_pz, pf.params->m_dx, pf.params->m_dy,pf.params->m_dz, pf.params->m_nu, dt);

        for(unsigned int i=0; i<pf.positions.size(); i++ )
        {
           pf.omega[i].x = tmp[i];
        }

        // y-dir
        for(unsigned int i=0; i<pf.positions.size(); i++ )
        {
           tmp[i] = pf.omega[i].y;
        }

        m_diffeqsolver3d->solve(tmp,pf.params->m_num_px,pf.params->m_num_py,pf.params->m_num_pz, pf.params->m_dx, pf.params->m_dy,pf.params->m_dz, pf.params->m_nu, dt);

        for(unsigned int i=0; i<pf.positions.size(); i++ )
        {
           pf.omega[i].y = tmp[i];
        }

        // z-dir
        for(unsigned int i=0; i<pf.positions.size(); i++ )
        {
           tmp[i] = pf.omega[i].z;
        }

        m_diffeqsolver3d->solve(tmp,pf.params->m_num_px,pf.params->m_num_py,pf.params->m_num_pz, pf.params->m_dx, pf.params->m_dy,pf.params->m_dz, pf.params->m_nu, dt);

        for(unsigned int i=0; i<pf.positions.size(); i++ )
        {
           pf.omega[i].z = tmp[i];
        }

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
