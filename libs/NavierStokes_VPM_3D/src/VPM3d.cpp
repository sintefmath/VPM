#include "VPM3d.hpp"

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

#define C_VPM 0.125

namespace VPM
{
    VPM3d::VPM3d(int argc, char** argv)
    {
        m_advection = std::make_shared<Split_Advection>();
        m_diffusion = std::make_shared<Split_Diffusion>();

        m_source = std::make_shared<Split_Source>();

    };
    VPM3d::~VPM3d()
    {
        m_advection.reset();
        m_diffusion.reset();
        m_source.reset();
    };

    bool VPM3d::new_stepsize(
            const double T,
            ParticleField & pf,
            double & delta_t
            )
    {
        bool ret = false;
        delta_t = C_VPM / std::max(pf.Linfty_omega, pf.Linfty_gradVelocity);
        if (pf.params->m_time+delta_t>T)
        {
            delta_t = T-pf.params->m_time;
            pf.params->m_time = T;
            ret=true;
        }
        else
        {
            pf.params->m_time += delta_t;
        }
        return ret;
    }

    void VPM3d::setStructure(std::shared_ptr<VPM::Structure> structure)
    {
        m_advection->setStructure(structure);
    }

    void VPM3d::run(
            const std::vector<Point3d> & positions,
            const std::vector<Point3d> & omega,
            const double T,
            const double dt,
            const std::shared_ptr<Parameters> params,
            const std::string outputFile,
            const int filename_count,
            const bool save_init
            )
    {


        // copy
        m_outputFile = outputFile;
        m_filename_count = filename_count;
        m_save_init = save_init;

        m_pf.omega = omega;
        m_pf.positions = positions;
        m_pf.regular_positions = positions;
        m_pf.params = params;
        m_pf.cartesianGrid = true;
        m_pf.velocity_correspondsTo_omega = false;

        m_pf.params->m_time = 0.;

        int stepcount = 0;
        double delta_t = 0.01;
        bool brexit = false;
        while (!brexit)
        {

// ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| //
            std::cerr<<"------------- Step " <<m_filename_count<<", time t = "<<m_pf.params->m_time<<" -----------\n";

// 1) Calculate step size from stability constraint
                std::cerr<<"  ---> calculate Linfty norm of the velocity field.\n";
                m_advection->getVelocityAndJacobian(m_pf,delta_t);
                std::cerr<<"    ---> max Linfty of velocity = "<<m_pf.Linfty_gradVelocity<<std::endl;
                if ( stepcount==0&&m_save_init )
                {
                    writeParticlesToFile(m_outputFile, m_filename_count, m_pf);
                    m_filename_count++;
                }

                if (stepcount==0)
                {
                    m_pf.Linfty_gradVelocity = C_VPM / delta_t;
                    m_pf.Linfty_omega = C_VPM / delta_t;
                }

                brexit = new_stepsize(T, m_pf, delta_t);

// 2) Diffusion
            m_diffusion->Euler_step( m_pf, delta_t);
// 3) Advection
            if (m_pf.params->m_order_ODEsolver==1)
            {
                m_advection->Euler_step( m_pf, delta_t );
            }
            else //if (m_pf.params->m_order_ODEsolver==2)
            {
                m_advection->RK2_step( m_pf, delta_t );
            }

            m_advection->getVelocityAndJacobian(m_pf,delta_t);
            writeParticlesToFile(m_outputFile, m_filename_count, m_pf);
            m_filename_count++;
            stepcount++;
            std::cerr<<"\n";

        }


    }
}
