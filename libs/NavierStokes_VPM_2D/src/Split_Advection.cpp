#include "Split_Advection.hpp"

#include "Utilities.hpp"
#include "PoissonSolver2D.hpp"
#include "Structure.hpp"

#include <assert.h>

namespace VPM
{

    Split_Advection::Split_Advection(double lambda)
        :
            m_lambda(lambda),
            m_structure_set(false)
    {
        PoissonSolver2D_BC bc = PoissonSolver2D_BC(
                POISSON2D_BC_NEUMANN,
                POISSON2D_BC_NEUMANN,
                POISSON2D_BC_NEUMANN,
                POISSON2D_BC_NEUMANN
                );
        m_poissonsolver2d = std::make_shared<ModulesHotel::PoissonSolver2D>(bc);
        m_redistribute = std::make_shared<Redistribute>();

    };

    Split_Advection::~Split_Advection()
    {
        m_poissonsolver2d.reset();
    };

    void Split_Advection::setStructure(std::shared_ptr<VPM::Structure> structure)
    {
        m_structure = structure;
        m_structure_set = true;
    }


    void Split_Advection::Euler_step(
            ParticleField & pf,
            const double dt
            )
    {

        if (!pf.cartesianGrid)
        {
            std::cerr<<"  ---> Split_Advection: Euler_step called on IRregular grid\n";
            exit(0);
            return;
        }

        std::cerr<<"  ---> Split_Advection: Euler_step\n";

        //p^{n+1} = p^n + dt*velocity(p^n);
        calculateVelocity(pf);

        for(unsigned int i=0; i<pf.positions.size(); i++ )
        {
           pf.positions[i] = pf.regular_positions[i] + dt*pf.velocity[i];
        }

        pf.cartesianGrid = false;
        pf.velocity_correspondsTo_omega = false;

        if (pf.params.m_remesh->m_isOn)
        {
            std::cerr<<"    --->Redistribute omega onto grid\n";
            m_redistribute->redistribute(pf);
        }

    }

    void Split_Advection::RK2_step(
            ParticleField & pf,
            const double dt
            )
    {

        if (!pf.cartesianGrid)
        {
            std::cerr<<"  ---> Split_Advection:: RK2_step called on IRregular grid\n";
            exit(0);
            return;
        }

        std::cerr<<"  ---> Split_Advection:: RK2_step\n";

        std::cerr<<"                          substep 1\n";
        //    p^* = p^n + dt/2*velocity(p^n)
        //p^{n+1} = p^n +   dt*velocity(p^*);
        calculateVelocity(pf);

        for(unsigned int i=0; i<pf.positions.size(); i++ )
        {
           pf.positions[i] = pf.regular_positions[i] + 0.5*dt*pf.velocity[i];
        }

        pf.cartesianGrid = false;
        pf.velocity_correspondsTo_omega = false;

        std::cerr<<"                          substep 2\n";
        calculateVelocity(pf);

        for(unsigned int i=0; i<pf.positions.size(); i++ )
        {
           pf.positions[i] = pf.regular_positions[i] + dt*pf.velocity[i];
        }

        pf.cartesianGrid = false;
        pf.velocity_correspondsTo_omega = false;

        if (pf.params.m_remesh->m_isOn)
        {
            std::cerr<<"    ---> Redistribute omega onto grid\n";
            m_redistribute->redistribute(pf);
        }

    }

    void Split_Advection::calculateVelocity(
            ParticleField & pf
            )
    {
        if (pf.velocity_correspondsTo_omega)
        {
            std::cerr<<"    ---> Velocity already calculated.\n";
            return;
        }

        // calculate potential
        std::vector<double> potential;
        if (pf.cartesianGrid)
        {
            std::cerr<<"    ---> Calculating velocity...";
            potential = m_poissonsolver2d->solve(pf.omega,pf.params.m_num_px,pf.params.m_num_py, pf.params.m_dx, pf.params.m_dy);
            std::cerr<<"done.\n";
        }
        else
        {
            std::cerr<<"    ---> Redistribute omega onto regular grid\n";
            std::vector<double> redistributed_omega;
            redistributed_omega.resize(pf.omega.size(), 0.);
            for (int i=0; i<pf.omega.size(); i++)
            {
                redistributed_omega[i] = pf.omega[i];
            }
            m_redistribute->redistribute(
                pf.positions,
                pf.regular_positions,
                pf.params,
                redistributed_omega,
                RED2D_G2P_BUILD
                );
            std::cerr<<"    ---> Calculating velocity...";
            potential = m_poissonsolver2d->solve(redistributed_omega,pf.params.m_num_px,pf.params.m_num_py, pf.params.m_dx, pf.params.m_dy);
            std::cerr<<"done.\n";
        }


        pf.velocity.clear();
        pf.velocity.resize(pf.omega.size(), Point2d(0.));

        for (int j=0; j<pf.params.m_regular_num_py; j++)
        {
            for (int i=0; i<pf.params.m_regular_num_px; i++)
            {
                // TODO boundary condition!!!
                // TODO boundary condition!!!
                // TODO boundary condition!!!
                int im = std::max(i-1,0);
                int ip = std::min(i+1,int(pf.params.m_regular_num_px-1));
                int jm = std::max(j-1,0);
                int jp = std::min(j+1,int(pf.params.m_regular_num_py-1));
                int i_j  = i +int(pf.params.m_regular_num_px)*j;
                int i_jm  = i +int(pf.params.m_regular_num_px)*jm;
                int i_jp  = i +int(pf.params.m_regular_num_px)*jp;
                int im_j  = im+int(pf.params.m_regular_num_px)*j;
                int ip_j  = ip+int(pf.params.m_regular_num_px)*j;
                pf.velocity[i_j] = pf.params.m_Uinfty + Point2d(
                                            (potential[i_jp] - potential[i_jm] )/(2.*pf.params.m_dy), // +dp/dy
                                           -(potential[ip_j] - potential[im_j] )/(2.*pf.params.m_dx)  // -dp/dx
                                        );
                if (m_structure_set)
                {
                    if (m_structure->isInside(pf.regular_positions[i_j], 0.))
                    {
                        //pf.velocity[i_j] = pf.velocity[i_j] / (1. + dt*m_lambda*1.);
                        //pf.velocity[i_j] = pf.velocity[i_j](1. - dt*m_lambda*1.);
                        pf.velocity[i_j] = Point2d(0,0);
                    }
                }
            }
        }

        // update vorticity
        if (m_structure_set)
        {
            for (int j=0; j<pf.params.m_regular_num_py; j++)
            {
                int jm = std::max(j-1,0);
                int jp = std::min(j+1,int(pf.params.m_num_py-1));
                for (int i=0; i<pf.params.m_regular_num_px; i++)
                {
                    int im = std::max(i-1,0);
                    int ip = std::min(i+1,int(pf.params.m_num_px-1));

                    int i_j   = i +int(pf.params.m_num_px)*j;

                    if (m_structure->isInside(pf.positions[i_j], 4*pf.params.m_dx))
                    {

                        int i_jm  = i +int(pf.params.m_num_px)*jm;
                        int i_jp  = i +int(pf.params.m_num_px)*jp;
                        int im_j  = im+int(pf.params.m_num_px)*j;
                        int ip_j  = ip+int(pf.params.m_num_px)*j;
//#ifdef LESS_DIFF
//                        pf.omega[i_j] += operator_curl(
//                                xi_i_jp*(SPEED_X-pf.velocity[i_jp].x),
//                                xi_i_jm*(SPEED_X-pf.velocity[i_jm].x),
//                                xi_ip_j*(SPEED_Y-pf.velocity[ip_j].y),
//                                xi_im_j*(SPEED_Y-pf.velocity[im_j].y),
//                                pf.params.m_dx,
//                                pf.params.m_dy
//                                );
//#else
//                        pf.omega[i_j] = operator_curl(
//                                (pf.velocity[i_jp].x+m_lambda*delta_t*xi_i_jp*SPEED_X)/(1.+m_lambda*delta_t*xi_i_jp),
//                                (pf.velocity[i_jm].x+m_lambda*delta_t*xi_i_jm*SPEED_X)/(1.+m_lambda*delta_t*xi_i_jm),
//                                (pf.velocity[ip_j].y+m_lambda*delta_t*xi_ip_j*SPEED_Y)/(1.+m_lambda*delta_t*xi_ip_j),
//                                (pf.velocity[im_j].y+m_lambda*delta_t*xi_im_j*SPEED_Y)/(1.+m_lambda*delta_t*xi_im_j),
//                                pf.params.m_dx,
//                                pf.params.m_dy
//                                );
//#endif
                        pf.omega[i_j] = operator_curl(
                                pf.velocity[i_jp].x,
                                pf.velocity[i_jm].x,
                                pf.velocity[ip_j].y,
                                pf.velocity[im_j].y,
                                pf.params.m_dx,
                                pf.params.m_dy
                                );
                    }
                }
            }
        }

        if (!pf.cartesianGrid)
        {
            if (m_structure_set)
            {
                std::cerr<<"    ---> Redistribute vorticity onto IRregular grid\n";
                m_redistribute->redistribute(
                    pf.regular_positions,
                    pf.positions,
                    pf.params,
                    pf.omega,
                    RED2D_G2P_USE
                            );
            }
            std::cerr<<"    ---> Redistribute velocity onto IRregular grid\n";
            m_redistribute->redistribute(
                pf.regular_positions,
                pf.positions,
                pf.params,
                pf.velocity,
                RED2D_G2P_USE
                );
        }

        pf.velocity_correspondsTo_omega = true;
        //for (int j=0; j<pf.params.m_regular_num_py; j++)
        //{
        //    for (int i=0; i<pf.params.m_regular_num_px; i++)
        //    {
        //        int i_j  = i +int(pf.params.m_regular_num_px)*j;
        //        std::cerr<<pf.velocity[i_j].x<<", ";
        //    }
        //    std::cerr<<"\n";
        //}
    }

    void Split_Advection::getLinftyGradVelocity(
            ParticleField & pf
            )
    {

        calculateVelocity( pf );

        pf.Linfty_gradVelocity = double(0);

        for (int j=1; j<pf.params.m_regular_num_py-1; j++)
        {
            for (int i=1; i<pf.params.m_regular_num_px-1; i++)
            {
                int im = /*std::max(*/i-1;//,0);
                int ip = /*std::min(*/i+1;//,int(pf.params.m_regular_num_px-1));
                int jm = /*std::max(*/j-1;//,0);
                int jp = /*std::min(*/j+1;//,int(pf.params.m_regular_num_py-1));
                int i_jm  = i +int(pf.params.m_regular_num_px)*jm;
                int i_jp  = i +int(pf.params.m_regular_num_px)*jp;
                int im_j  = im+int(pf.params.m_regular_num_px)*j;
                int ip_j  = ip+int(pf.params.m_regular_num_px)*j;

                pf.Linfty_gradVelocity = std::max(pf.Linfty_gradVelocity, std::abs(
                            (pf.velocity[ip_j].x - pf.velocity[im_j].x )
                            /(2.*pf.params.m_dx)) );
                pf.Linfty_gradVelocity = std::max(pf.Linfty_gradVelocity, std::abs(
                            (pf.velocity[i_jp].x - pf.velocity[i_jm].x )
                            /(2.*pf.params.m_dy)) );
                pf.Linfty_gradVelocity = std::max(pf.Linfty_gradVelocity, std::abs(
                            (pf.velocity[ip_j].y - pf.velocity[im_j].y )
                            /(2.*pf.params.m_dx)) );
                pf.Linfty_gradVelocity = std::max(pf.Linfty_gradVelocity, std::abs(
                            (pf.velocity[i_jp].y - pf.velocity[i_jm].y )
                            /(2.*pf.params.m_dy)) );

            }
        }

    }

}
