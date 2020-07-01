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

#define C_VPM 0.125

namespace VPM
{
    VPM2d::VPM2d(int argc, char** argv)
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

    void VPM2d::setStructure(std::shared_ptr<VPM::Structure> structure)
    {
        m_advection->setStructure(structure);
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
            m_filename_count++;
            stepcount++;
            writeParticlesToFile(m_outputFile, m_filename_count, m_pf);
            if (onestep) {
                brexit = true;
            }

        }


    }
}
