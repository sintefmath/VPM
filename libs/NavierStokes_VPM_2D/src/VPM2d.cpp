#include "VPM2d.hpp"

#include "Split_Advection.hpp"
#include "Split_Diffusion.hpp"
#include "Split_Source.hpp"

#include "Utilities.hpp"

#include <iostream>
#include <fstream>
#include <math.h>

#include <algorithm>

#include <time.h>

#include <petscksp.h>

#define C_VPM 0.75

namespace VPM
{
    VPM2d::VPM2d(int argc, char** argv)
        :
        m_structure_set(false)
    {
        m_advection = std::make_shared<Split_Advection>();
        m_diffusion = std::make_shared<Split_Diffusion>();

        m_source = std::make_shared<Split_Source>();

    };
    VPM2d::~VPM2d()
    {
        m_advection.reset();
        m_diffusion.reset();
        m_source.reset();
    };

    bool VPM2d::new_stepsize(
            const double final_time,
            ParticleField & pf,
            double & delta_t
            )
    {
        bool ret = false;
        delta_t = C_VPM / pf.Linfty_gradVelocity;
        if (pf.time+delta_t>final_time)
        {
            delta_t = final_time-pf.time;
            pf.time = final_time;
            ret=true;
        }
        else
        {
            pf.time += delta_t;
        }
        return ret;
    }

    // implements equation (32) in https://www-ljk.imag.fr/membres/Chloe.Mimeau/paper_IJNMF.pdf
    Point2d VPM2d::getForces(
            const ParticleField & pf,
            const ParticleField & pf_old,
            const double & delta_t
            )
    {
        Point2d F = Point2d(0,0);
        Point2d origo;
        m_structure->getOrigo(origo);
        double Reynolds_inv = pf.params.m_nu/(m_structure->getCharacteristicLength()*L2Norm(pf.params.m_Uinfty));
        if (m_structure_set)
        {
            if (!pf.cartesianGrid)
            {
                std::cerr<<"  ---> calculating forces called on IRregular grid\n";
                exit(0);
                return Point2d();
            }
#define forcesversion2
#ifdef forcesversion2
            int c = 3;
            for (int j=c; j<pf.params.m_regular_num_py-c; j++)
            {
                for (int i=c; i<pf.params.m_regular_num_px-c; i++)
                {
                    int i_j  = i +int(pf.params.m_regular_num_px)*j;
                    //if (!m_structure->isInside(pf.regular_positions[i_j], 0.) && m_structure->isInside(pf.regular_positions[i_j], .4))
                    //{
                        F -= (pf.velocity[i_j]-pf_old.velocity[i_j])/delta_t*pf.params.m_vol;
                    //}
                }
            }
            for (int i=c-1; i<pf.params.m_regular_num_px-c+1; i++)
            {
                for (int ul = 0; ul<2; ul++)
                {
                    double sign=+1;
                    int j=c-1;
                    if (ul>0)
                    {
                        j=pf.params.m_regular_num_py-c;
                        sign=-1;
                    }
                    int i_j   = i  +int(pf.params.m_regular_num_px)*j;
                    int i_jm  = i  +int(pf.params.m_regular_num_px)*(j-1);
                    int i_jp  = i  +int(pf.params.m_regular_num_px)*(j+1);
                    int ip_jp = i+1+int(pf.params.m_regular_num_px)*(j+1);
                    int ip_jm = i+1+int(pf.params.m_regular_num_px)*(j-1);
                    int im_jp = i-1+int(pf.params.m_regular_num_px)*(j+1);
                    int im_jm = i-1+int(pf.params.m_regular_num_px)*(j-1);
                    int im_j  = i-1+int(pf.params.m_regular_num_px)*j;
                    int ip_j  = i+1+int(pf.params.m_regular_num_px)*j;
                    double x = (pf.regular_positions[i_j].x-origo.x);
                    double y = (pf.regular_positions[i_j].y-origo.y);
                    double uy = (pf.velocity[i_jp].x - pf.velocity[i_jm].x)/(2*pf.params.m_dy);
                    double vx = (pf.velocity[ip_j].y - pf.velocity[im_j].y)/(2*pf.params.m_dx);
                    double vy = (pf.velocity[i_jp].y - pf.velocity[i_jm].y)/(2*pf.params.m_dy);
                    double uxx = (pf.velocity[ip_j].x - 2*pf.velocity[i_j].x + pf.velocity[im_j].x)/pf.params.m_vol;
                    double uyy = (pf.velocity[i_jp].x - 2*pf.velocity[i_j].x + pf.velocity[i_jm].x)/pf.params.m_vol;
                    double vxy = (pf.velocity[ip_jp].y - pf.velocity[ip_jm].y - pf.velocity[im_jp].y + pf.velocity[im_jm].y)/pf.params.m_vol;
                    double ut = (pf.velocity[i_j].x-pf_old.velocity[i_j].x)/delta_t;
                    F.x += sign*pf.params.m_dx*(
                           pf.velocity[i_j].x*pf.velocity[i_j].y
                         + pf.velocity[i_j].y*pf.omega[i_j]*y
                         - y*ut
                         + Reynolds_inv*(2*uxx + uyy + vxy)*y
                         - Reynolds_inv*(uy + vx));
                    F.y += sign*pf.params.m_dx*(
                           0.5*(pf.velocity[i_j].y*pf.velocity[i_j].y - pf.velocity[i_j].x*pf.velocity[i_j].x)
                         - pf.velocity[i_j].y*pf.omega[i_j]*x
                         + x*ut
                         - Reynolds_inv*(2*uxx + uyy + vxy)*x
                         - 2*Reynolds_inv*vy);
                }
            }
            for (int j=c-1; j<pf.params.m_regular_num_py-c+1; j++)
            {
                for (int ul = 0; ul<2; ul++)
                {
                    double sign=-1;
                    int i=c-1;
                    if (ul>0)
                    {
                        i=pf.params.m_regular_num_px-c;
                        sign=+1;
                    }
                    int i_j   = i  +int(pf.params.m_regular_num_px)*j;
                    int i_jm  = i  +int(pf.params.m_regular_num_px)*(j-1);
                    int i_jp  = i  +int(pf.params.m_regular_num_px)*(j+1);
                    int ip_jp = i+1+int(pf.params.m_regular_num_px)*(j+1);
                    int ip_jm = i+1+int(pf.params.m_regular_num_px)*(j-1);
                    int im_jp = i-1+int(pf.params.m_regular_num_px)*(j+1);
                    int im_jm = i-1+int(pf.params.m_regular_num_px)*(j-1);
                    int im_j  = i-1+int(pf.params.m_regular_num_px)*j;
                    int ip_j  = i+1+int(pf.params.m_regular_num_px)*j;
                    double x = (pf.regular_positions[i_j].x-origo.x);
                    double y = (pf.regular_positions[i_j].y-origo.y);
                    double ux = (pf.velocity[ip_j].x - pf.velocity[im_j].x)/(2*pf.params.m_dx);
                    double uy = (pf.velocity[i_jp].x - pf.velocity[i_jm].x)/(2*pf.params.m_dy);
                    double vx = (pf.velocity[ip_j].y - pf.velocity[im_j].y)/(2*pf.params.m_dx);
                    double vxx = (pf.velocity[ip_j].y - 2*pf.velocity[i_j].y + pf.velocity[im_j].y)/pf.params.m_vol;
                    double vyy = (pf.velocity[i_jp].y - 2*pf.velocity[i_j].y + pf.velocity[i_jm].y)/pf.params.m_vol;
                    double uxy = (pf.velocity[ip_jp].x - pf.velocity[ip_jm].x - pf.velocity[im_jp].x + pf.velocity[im_jm].x)/pf.params.m_vol;
                    double vt = (pf.velocity[i_j].y-pf_old.velocity[i_j].y)/delta_t;
                    F.x += sign*pf.params.m_dx*(
                           0.5*(pf.velocity[i_j].y*pf.velocity[i_j].y-pf.velocity[i_j].x*pf.velocity[i_j].x)
                         - pf.velocity[i_j].x*pf.omega[i_j]*y
                         - y*vt
                         + Reynolds_inv*(2*vyy + vxx + uxy)*y
                         + 2*Reynolds_inv*ux);
                    F.y += sign*pf.params.m_dx*(
                         - pf.velocity[i_j].x*pf.velocity[i_j].y
                         + pf.velocity[i_j].x*pf.omega[i_j]*x
                         + x*vt
                         - Reynolds_inv*(2*vyy + vxx + uxy)*x
                         + Reynolds_inv*(vx + uy));
                }
            }
#else
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
                    int i_j   = i +int(pf.params.m_regular_num_px)*j;
                    int i_jm  = i +int(pf.params.m_regular_num_px)*jm;
                    int i_jp  = i +int(pf.params.m_regular_num_px)*jp;
                    int im_j  = im+int(pf.params.m_regular_num_px)*j;
                    int ip_j  = ip+int(pf.params.m_regular_num_px)*j;

                    double nabla2omega = pf.omega[ip_j] + pf.omega[im_j] -4.*pf.omega[i_j] + pf.omega[i_jp] + pf.omega[i_jm];
                    // The following line is correct, but it canceles with the line below
                    // nabla2omega /= pf.params.m_vol;

                    Fx -=  pf.params.m_nu*nabla2omega*(pf.regular_positions[i_j].y-origo.y);
                    Fy +=  pf.params.m_nu*nabla2omega*(pf.regular_positions[i_j].x-origo.x);
                    // Fx *= pf.params.m_vol;
                    // Fy *= pf.params.m_vol;
                }
            }
#endif
        }

        return F;
    }

    void VPM2d::setStructure(std::shared_ptr<VPM::Structure> structure)
    {
        m_advection->setStructure(structure);
        m_structure = structure;
        m_structure_set = true;
    }

    void VPM2d::run(
            ParticleField & pf,
            //const std::vector<Point2d> & positions,
            //const std::vector<double> & omega,
            const double final_time,
            //const std::shared_ptr<Parameters> params,
            const std::string outputFile,
            const int filename_count,
            const bool save_init,
            const bool onestep
            )
    {

        m_pf = pf;


        // copy
        m_outputFile = outputFile;
        m_filename_count = filename_count;
        m_save_init = save_init;

        std::vector<double> ft, fx, fy;

        int stepcount = 0;
        double delta_t = 0.01;
        bool brexit = false;
        while (!brexit)
        {

// ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| //
            std::cerr<<"------------- Step " <<m_filename_count<<", time t = "<<m_pf.time<<" -----------\n";
// 1) Calculate step size from stability constraint
            std::cerr<<"  ---> calculate Linfty norm of the velocity field.\n";
            m_advection->getLinftyGradVelocity(m_pf);
            std::cerr<<"    ---> max Linfty of velocity = "<<m_pf.Linfty_gradVelocity<<std::endl;
            if ( stepcount==0&&m_save_init )
            {
                writeParticlesToFile(m_outputFile, m_filename_count, m_pf);
                m_filename_count++;
            }

            //if (stepcount==0)
            //{
            //    m_pf.Linfty_gradVelocity = C_VPM / delta_t;
            //}

            brexit = new_stepsize(final_time, m_pf, delta_t);

            m_pf_old = pf;

// 2) Diffusion
            m_diffusion->Euler_step( m_pf, delta_t);
// 3) Advection
            if (m_pf.params.m_order_ODEsolver==1)
            {
                m_advection->Euler_step( m_pf, delta_t );
            }
            else //if (m_pf.params.m_order_ODEsolver==2)
            {
                m_advection->RK2_step( m_pf, delta_t );
            }

            m_advection->calculateVelocity(m_pf);
            Point2d F = getForces(m_pf, m_pf_old, delta_t);

            std::cerr<<"Forces="<<F.x<<", "<<F.y<<std::endl;

            ft.push_back(m_pf.time);
            fx.push_back(F.x);
            fy.push_back(F.y);


            m_filename_count++;
            stepcount++;
            if(stepcount%100==0)
            {
                double ForceToCoeff = 2./m_pf.params.m_Uinfty.x/m_pf.params.m_Uinfty.x/m_structure->getCharacteristicLength();
                writeParticlesToFile(m_outputFile, m_filename_count, m_pf);

                std::string t = std::to_string(m_filename_count);
                std::ofstream myfile_ft;
                std::ofstream myfile_fx;
                std::ofstream myfile_fy;
                std::ofstream myfile_dx;
                std::ofstream myfile_dy;
                myfile_ft.open(m_outputFile+"_t"+t+"_ft.dat", std::ofstream::binary);
                myfile_fx.open(m_outputFile+"_t"+t+"_fx.dat", std::ofstream::binary);
                myfile_fy.open(m_outputFile+"_t"+t+"_fy.dat", std::ofstream::binary);
                myfile_dx.open(m_outputFile+"_t"+t+"_dx.dat", std::ofstream::binary);
                myfile_dy.open(m_outputFile+"_t"+t+"_dy.dat", std::ofstream::binary);
                double tmp;
                for (unsigned int i=0; i<fx.size(); i++)
                {
#ifndef forcesversion1
                    myfile_ft.write(reinterpret_cast<const char*>(&ft[i]), sizeof(double));
                    tmp = fx[i];
                    myfile_fx.write(reinterpret_cast<const char*>(&tmp), sizeof(double));
                    tmp *= ForceToCoeff;
                    myfile_dx.write(reinterpret_cast<const char*>(&tmp), sizeof(double));
                    tmp = fy[i];
                    myfile_fy.write(reinterpret_cast<const char*>(&tmp), sizeof(double));
                    tmp *= ForceToCoeff;
                    myfile_dy.write(reinterpret_cast<const char*>(&tmp), sizeof(double));
#else
                    myfile_ft.write(reinterpret_cast<const char*>(&ft[i]), sizeof(double));
                    myfile_fx.write(reinterpret_cast<const char*>(&fx[i]), sizeof(double));
                    myfile_fy.write(reinterpret_cast<const char*>(&fy[i]), sizeof(double));
                    myfile_dx.write(reinterpret_cast<const char*>(&dx[i]), sizeof(double));
                    myfile_dy.write(reinterpret_cast<const char*>(&dy[i]), sizeof(double));
#endif
                }
                myfile_ft.close();
                myfile_fx.close();
                myfile_fy.close();
                myfile_dx.close();
                myfile_dy.close();

            }
            if (onestep) {
                brexit = true;
            }

        }


    }
}
