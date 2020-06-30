#include "Split_Advection.hpp"

#include "Utilities.hpp"
#include "PoissonSolver3D.hpp"
#include "FlattenArrayIndex.hpp"
#include "Structure.hpp"

#include <assert.h>

#define TRANSPOSE_SCHEME

namespace VPM
{

    Split_Advection::Split_Advection(double lambda)
        :
            m_lambda(lambda),
            m_structure_set(false)
    {
        PoissonSolver3D_BC bc = PoissonSolver3D_BC(
                POISSON3D_BC_NEUMANN,
                POISSON3D_BC_NEUMANN,
                POISSON3D_BC_NEUMANN,
                POISSON3D_BC_NEUMANN,
                POISSON3D_BC_NEUMANN,
                POISSON3D_BC_NEUMANN
                );
        m_poissonsolver3d = std::make_shared<ModulesHotel::PoissonSolver3D>(bc);
        m_redistribute = std::make_shared<Redistribute>();

    };

    Split_Advection::~Split_Advection()
    {
        m_poissonsolver3d.reset();
        m_redistribute.reset();
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
        getVelocityAndJacobian(pf,dt);

        for(unsigned int i=0; i<pf.positions.size(); i++ )
        {
           pf.positions[i] = pf.regular_positions[i] + dt*pf.velocity[i];
           pf.omega[i] = pf.omega[i] + dt*pf.Jvel_x_omega[i];
        }

        pf.cartesianGrid = false;
        pf.velocity_correspondsTo_omega = false;

        if (pf.params->m_remesh->m_isOn)
        {
            std::cerr<<"    ---> Redistribute omega onto grid\n";
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
            std::cerr<<"  ---> Split_Advection: RK2_step called on IRregular grid\n";
            exit(0);
            return;
        }

        std::cerr<<"  ---> Split_Advection: RK2_step\n";

        std::cerr<<"                          substep 1\n";
        //    p^* = p^n + dt/2*velocity(p^n)
        //p^{n+1} = p^n +   dt*velocity(p^*);
        getVelocityAndJacobian(pf,dt);

        std::vector<Point3d> omega_tmp;
        omega_tmp.clear();
        omega_tmp.resize(pf.positions.size(), Point3d(0.));

        for(unsigned int i=0; i<pf.positions.size(); i++ )
        {
           pf.positions[i] = pf.regular_positions[i] + 0.5*dt*pf.velocity[i];
           omega_tmp[i] = pf.omega[i];
           pf.omega[i] = pf.omega[i] + 0.5*dt*pf.Jvel_x_omega[i];
        }

        pf.cartesianGrid = false;
        pf.velocity_correspondsTo_omega = false;

        std::cerr<<"                          substep 2\n";
        getVelocityAndJacobian(pf,dt);

        for(unsigned int i=0; i<pf.positions.size(); i++ )
        {
           pf.positions[i] = pf.regular_positions[i] + dt*pf.velocity[i];
           pf.omega[i] = omega_tmp[i] + dt*pf.Jvel_x_omega[i];
        }

        pf.cartesianGrid = false;
        pf.velocity_correspondsTo_omega = false;

        if (pf.params->m_remesh->m_isOn)
        {
            std::cerr<<"    ---> Redistribute omega onto grid\n";
            m_redistribute->redistribute(pf);
        }

    }

    void Split_Advection::getVelocityAndJacobian(
            ParticleField & pf,
            const double dt
            )
    {
        if (pf.velocity_correspondsTo_omega)
        {
            std::cerr<<"    ---> Velocity already calculated.\n";
            return;
        }

        // calculate potential
        std::vector<double> potentialx;
        std::vector<double> potentialy;
        std::vector<double> potentialz;

        if (pf.cartesianGrid)
        {
            std::cerr<<"    ---> Calculating velocity x-dir...";
            std::vector<double> w;
            w.resize(pf.params->m_N);
            for (int i=0; i<pf.params->m_N; i++)
            {
                w[i] = pf.omega[i].x;
            }
            potentialx = m_poissonsolver3d->solve(w,pf.params->m_num_px,pf.params->m_num_py,pf.params->m_num_pz, pf.params->m_dx, pf.params->m_dy, pf.params->m_dz);
            std::cerr<<"done.\n";

            std::cerr<<"    ---> Calculating velocity y-dir...";
            for (int i=0; i<pf.params->m_N; i++)
            {
                w[i] = pf.omega[i].y;
            }
            potentialy = m_poissonsolver3d->solve(w,pf.params->m_num_px,pf.params->m_num_py,pf.params->m_num_pz, pf.params->m_dx, pf.params->m_dy, pf.params->m_dz);
            std::cerr<<"done.\n";

            std::cerr<<"    ---> Calculating velocity z-dir...";
            for (int i=0; i<pf.params->m_N; i++)
            {
                w[i] = pf.omega[i].z;
            }
            potentialz = m_poissonsolver3d->solve(w,pf.params->m_num_px,pf.params->m_num_py,pf.params->m_num_pz, pf.params->m_dx, pf.params->m_dy, pf.params->m_dz);
            std::cerr<<"done.\n";
        }
        else
        {
            std::cerr<<"    ---> Redistribute omega onto regular grid\n";
            std::vector<Point3d> redistributed_omega;
            redistributed_omega.resize(pf.omega.size());
            for (int i=0; i<pf.omega.size(); i++)
            {
                redistributed_omega[i] = pf.omega[i];
            }
            m_redistribute->redistribute(
                pf.positions,
                pf.regular_positions,
                pf.params,
                redistributed_omega,
                RED3D_G2P_BUILD
                );

            std::vector<double> tmp;
            tmp.resize(redistributed_omega.size());

            for (int i=0; i<redistributed_omega.size(); i++)
            {
                tmp[i] = redistributed_omega[i].x;
            }
            std::cerr<<"    ---> Calculating velocity x-dir...";
            potentialx = m_poissonsolver3d->solve(tmp,pf.params->m_num_px,pf.params->m_num_py,pf.params->m_num_pz, pf.params->m_dx, pf.params->m_dy, pf.params->m_dz);
            std::cerr<<"done.\n";

            for (int i=0; i<pf.omega.size(); i++)
            {
                tmp[i] = redistributed_omega[i].y;
            }
            std::cerr<<"    ---> Calculating velocity y-dir...";
            potentialy = m_poissonsolver3d->solve(tmp,pf.params->m_num_px,pf.params->m_num_py,pf.params->m_num_pz, pf.params->m_dx, pf.params->m_dy, pf.params->m_dz);
            std::cerr<<"done.\n";

            for (int i=0; i<pf.omega.size(); i++)
            {
                tmp[i] = redistributed_omega[i].z;
            }
            std::cerr<<"    ---> Calculating velocity z-dir...";
            potentialz = m_poissonsolver3d->solve(tmp,pf.params->m_num_px,pf.params->m_num_py,pf.params->m_num_pz, pf.params->m_dx, pf.params->m_dy, pf.params->m_dz);
            std::cerr<<"done.\n";
        }


        // create vector of correct size the first time
        if (pf.velocity.size()!=pf.omega.size())
        {
            pf.velocity.clear();
            pf.velocity.resize(pf.omega.size(), Point3d(0.));
        }

        for (int k=0; k<pf.params->m_regular_num_pz; k++)
        {
            int km = std::max(k-1,0);
            int kp = std::min(k+1,int(pf.params->m_regular_num_pz-1));
            for (int j=0; j<pf.params->m_regular_num_py; j++)
            {
                int jm = std::max(j-1,0);
                int jp = std::min(j+1,int(pf.params->m_regular_num_py-1));
                for (int i=0; i<pf.params->m_regular_num_px; i++)
                {
                    // TODO boundary condition!!!
                    // TODO boundary condition!!!
                    // TODO boundary condition!!!
                    int im = std::max(i-1,0);
                    int ip = std::min(i+1,int(pf.params->m_regular_num_px-1));

                    int i_j_k  = flat(i,j,k,int(pf.params->m_regular_num_px),int(pf.params->m_regular_num_py));

                    int im_j_k  = flat(im,j,k,int(pf.params->m_regular_num_px),int(pf.params->m_regular_num_py));
                    int ip_j_k  = flat(ip,j,k,int(pf.params->m_regular_num_px),int(pf.params->m_regular_num_py));
                    int i_jm_k  = flat(i,jm,k,int(pf.params->m_regular_num_px),int(pf.params->m_regular_num_py));
                    int i_jp_k  = flat(i,jp,k,int(pf.params->m_regular_num_px),int(pf.params->m_regular_num_py));
                    int i_j_km  = flat(i,j,km,int(pf.params->m_regular_num_px),int(pf.params->m_regular_num_py));
                    int i_j_kp  = flat(i,j,kp,int(pf.params->m_regular_num_px),int(pf.params->m_regular_num_py));

                    pf.velocity[i_j_k] = pf.params->m_Uinfty + operator_curl(
                        potentialx[i_j_km], potentialx[i_j_kp],
                        potentialx[i_jm_k], potentialx[i_jp_k],
                        potentialy[im_j_k], potentialy[ip_j_k],
                        potentialy[i_j_km], potentialy[i_j_kp],
                        potentialz[im_j_k], potentialz[ip_j_k],
                        potentialz[i_jm_k], potentialz[i_jp_k],
                        pf.params->m_dx, pf.params->m_dy, pf.params->m_dz
                        );
                    if (m_structure_set)
                    {
                        if (m_structure->isInside(pf.regular_positions[i_j_k], 0.))
                        {
                            //pf.velocity[i_j_k] = pf.velocity[i_j_k] / (1. + dt*m_lambda*1.);
                            //pf.velocity[i_j_k] = pf.velocity[i_j_k](1. - dt*m_lambda*1.);
                            pf.velocity[i_j_k] = Point3d(0,0,0);
                        }
                    }
                }
            }
        }

        // update vorticity
        if (m_structure_set)
        {
            for (int k=0; k<pf.params->m_regular_num_pz; k++)
            {
                int km = std::max(k-1,0);
                int kp = std::min(k+1,int(pf.params->m_num_pz-1));
                for (int j=0; j<pf.params->m_regular_num_py; j++)
                {
                    int jm = std::max(j-1,0);
                    int jp = std::min(j+1,int(pf.params->m_num_py-1));
                    for (int i=0; i<pf.params->m_regular_num_px; i++)
                    {
                        int im = std::max(i-1,0);
                        int ip = std::min(i+1,int(pf.params->m_num_px-1));

                        int i_j_k = flat(i,j,k, int(pf.params->m_num_px), int(pf.params->m_num_py));

                        if (m_structure->isInside(pf.positions[i_j_k], 4*pf.params->m_dx))
                        {

                            int im_j_k = flat(im,j,k, int(pf.params->m_num_px), int(pf.params->m_num_py));
                            int ip_j_k = flat(ip,j,k, int(pf.params->m_num_px), int(pf.params->m_num_py));
                            int i_jm_k = flat(i,jm,k, int(pf.params->m_num_px), int(pf.params->m_num_py));
                            int i_jp_k = flat(i,jp,k, int(pf.params->m_num_px), int(pf.params->m_num_py));
                            int i_j_km = flat(i,j,km, int(pf.params->m_num_px), int(pf.params->m_num_py));
                            int i_j_kp = flat(i,j,kp, int(pf.params->m_num_px), int(pf.params->m_num_py));

//#ifdef LESS_DIFF
//                    pf.omega[i_j] += operator_curl(
//                            -xi_i_jp*pf.velocity[i_jp].x,
//                            -xi_i_jm*pf.velocity[i_jm].x,
//                            -xi_ip_j*pf.velocity[ip_j].y,
//                            -xi_im_j*pf.velocity[im_j].y,
//                            pf.params->m_dx,
//                            pf.params->m_dy
//                            );
//#else
//                    pf.omega[i_j] = operator_curl(
//                            pf.velocity[i_jp].x/(1.+m_lambda*delta_t*xi_i_jp),
//                            pf.velocity[i_jm].x/(1.+m_lambda*delta_t*xi_i_jm),
//                            pf.velocity[ip_j].y/(1.+m_lambda*delta_t*xi_ip_j),
//                            pf.velocity[im_j].y/(1.+m_lambda*delta_t*xi_im_j),
//                            pf.params->m_dx,
//                            pf.params->m_dy
//                            );
//#endif

                            pf.omega[i_j_k] = operator_curl(
                                    pf.velocity[i_j_km].x, pf.velocity[i_j_kp].x,
                                    pf.velocity[i_jm_k].x, pf.velocity[i_jp_k].x,
                                    pf.velocity[im_j_k].y, pf.velocity[ip_j_k].y,
                                    pf.velocity[i_j_km].y, pf.velocity[i_j_kp].y,
                                    pf.velocity[im_j_k].z, pf.velocity[ip_j_k].z,
                                    pf.velocity[i_jm_k].z, pf.velocity[i_jp_k].z,
                                    pf.params->m_dx, pf.params->m_dy, pf.params->m_dz
                                    );
                        }
                    }
                }
            }
        }

        calcLinftyGradVelocity(pf);
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
                    RED3D_G2P_USE
                            );
            }
            std::cerr<<"    ---> Redistribute velocity onto IRregular grid\n";
            m_redistribute->redistribute(
                pf.regular_positions,
                pf.positions,
                pf.params,
                pf.velocity,
                RED3D_G2P_USE
                        );
            std::cerr<<"    ---> Redistribute Dvel*omega onto IRregular grid\n";
            m_redistribute->redistribute(
                pf.regular_positions,
                pf.positions,
                pf.params,
                pf.Jvel_x_omega,
                RED3D_G2P_USE
                        );
                // TODO: in 3D project to be divergence free ?!
                // TODO: in 3D project to be divergence free ?!
                // TODO: in 3D project to be divergence free ?!
        }

        pf.velocity_correspondsTo_omega = true;

    }

    void Split_Advection::calcLinftyGradVelocity(
            ParticleField & pf
            )
    {

        pf.Linfty_gradVelocity = double(0);
        pf.Linfty_omega = double(0);
        if (pf.Jvel_x_omega.size()!=pf.omega.size())
        {
            pf.Jvel_x_omega.clear();
            pf.Jvel_x_omega.resize(pf.omega.size(), Point3d(0.));
        }

        for (int k=0; k<pf.params->m_regular_num_pz; k++)
        {
            int km = std::max(k-1,0);
            int kp = std::min(k+1,int(pf.params->m_regular_num_pz-1));
            for (int j=0; j<pf.params->m_regular_num_py; j++)
            {
                int jm = std::max(j-1,0);
                int jp = std::min(j+1,int(pf.params->m_regular_num_py-1));
                for (int i=0; i<pf.params->m_regular_num_px; i++)
                {
                    int im = std::max(i-1,0);
                    int ip = std::min(i+1,int(pf.params->m_regular_num_px-1));

                    int i_j_k  = flat(i,j,k,int(pf.params->m_regular_num_px),int(pf.params->m_regular_num_py));

                    int im_j_k  = flat(im,j,k,int(pf.params->m_regular_num_px),int(pf.params->m_regular_num_py));
                    int ip_j_k  = flat(ip,j,k,int(pf.params->m_regular_num_px),int(pf.params->m_regular_num_py));
                    int i_jm_k  = flat(i,jm,k,int(pf.params->m_regular_num_px),int(pf.params->m_regular_num_py));
                    int i_jp_k  = flat(i,jp,k,int(pf.params->m_regular_num_px),int(pf.params->m_regular_num_py));
                    int i_j_km  = flat(i,j,km,int(pf.params->m_regular_num_px),int(pf.params->m_regular_num_py));
                    int i_j_kp  = flat(i,j,kp,int(pf.params->m_regular_num_px),int(pf.params->m_regular_num_py));

                    double velx_dx = (pf.velocity[ip_j_k].x - pf.velocity[im_j_k].x )/(2.*pf.params->m_dx);
                    pf.Linfty_gradVelocity = std::max(pf.Linfty_gradVelocity, std::abs(velx_dx));
                    double velx_dy = (pf.velocity[i_jp_k].x - pf.velocity[i_jm_k].x )/(2.*pf.params->m_dy);
                    pf.Linfty_gradVelocity = std::max(pf.Linfty_gradVelocity, std::abs(velx_dy));
                    double velx_dz = (pf.velocity[i_j_kp].x - pf.velocity[i_j_km].x )/(2.*pf.params->m_dz);
                    pf.Linfty_gradVelocity = std::max(pf.Linfty_gradVelocity, std::abs(velx_dz));

                    double vely_dx = (pf.velocity[ip_j_k].y - pf.velocity[im_j_k].y )/(2.*pf.params->m_dx);
                    pf.Linfty_gradVelocity = std::max(pf.Linfty_gradVelocity, std::abs(vely_dx));
                    double vely_dy = (pf.velocity[i_jp_k].y - pf.velocity[i_jm_k].y )/(2.*pf.params->m_dy);
                    pf.Linfty_gradVelocity = std::max(pf.Linfty_gradVelocity, std::abs(vely_dy));
                    double vely_dz = (pf.velocity[i_j_kp].y - pf.velocity[i_j_km].y )/(2.*pf.params->m_dz);
                    pf.Linfty_gradVelocity = std::max(pf.Linfty_gradVelocity, std::abs(vely_dz));
                                

                    double velz_dx = (pf.velocity[ip_j_k].z - pf.velocity[im_j_k].z )/(2.*pf.params->m_dx);
                    pf.Linfty_gradVelocity = std::max(pf.Linfty_gradVelocity, std::abs(velz_dx));
                    double velz_dy = (pf.velocity[i_jp_k].z - pf.velocity[i_jm_k].z )/(2.*pf.params->m_dy);
                    pf.Linfty_gradVelocity = std::max(pf.Linfty_gradVelocity, std::abs(velz_dy));
                    double velz_dz = (pf.velocity[i_j_kp].z - pf.velocity[i_j_km].z )/(2.*pf.params->m_dz);
                    pf.Linfty_gradVelocity = std::max(pf.Linfty_gradVelocity, std::abs(velz_dz));

                    pf.Linfty_omega = std::max(pf.Linfty_omega,std::abs(pf.omega[i_j_k].x));
                    pf.Linfty_omega = std::max(pf.Linfty_omega,std::abs(pf.omega[i_j_k].y));
                    pf.Linfty_omega = std::max(pf.Linfty_omega,std::abs(pf.omega[i_j_k].z));

#ifdef TRANSPOSE_SCHEME
                    pf.Jvel_x_omega[i_j_k].x = velx_dx*pf.omega[i_j_k].x
                                             + vely_dx*pf.omega[i_j_k].y
                                             + velz_dx*pf.omega[i_j_k].z;
                    pf.Jvel_x_omega[i_j_k].y = velx_dy*pf.omega[i_j_k].x
                                             + vely_dy*pf.omega[i_j_k].y
                                             + velz_dy*pf.omega[i_j_k].z;
                    pf.Jvel_x_omega[i_j_k].z = velx_dz*pf.omega[i_j_k].x
                                             + vely_dz*pf.omega[i_j_k].y
                                             + velz_dz*pf.omega[i_j_k].z;
#else
                    pf.Jvel_x_omega[i_j_k].x = velx_dx*pf.omega[i_j_k].x
                                             + velx_dy*pf.omega[i_j_k].y
                                             + velx_dz*pf.omega[i_j_k].z;
                    pf.Jvel_x_omega[i_j_k].y = vely_dx*pf.omega[i_j_k].x
                                             + vely_dy*pf.omega[i_j_k].y
                                             + vely_dz*pf.omega[i_j_k].z;
                    pf.Jvel_x_omega[i_j_k].z = velz_dx*pf.omega[i_j_k].x
                                             + velz_dy*pf.omega[i_j_k].y
                                             + velz_dz*pf.omega[i_j_k].z;
#endif

                }
            }
        }

    }

}
