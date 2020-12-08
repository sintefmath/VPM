#pragma once
#include "Particles2d.hpp"
#include "Parameters.hpp"
#include "InitialConditions.hpp"
#include "Structure.hpp"

#include <vector>
#include <memory>

namespace VPM
{
    class Split_Advection;
    class Split_Diffusion;
    class Split_Source;
}

namespace VPM
{
    class VPM2d
    {
    public:
        VPM2d(int argc, char** argv);
        ~VPM2d();

        void run(
            ParticleField & pf,
            //const std::vector<Point2d> & positions,
            //const std::vector<double> & omega,
            const double final_time,
            //const std::shared_ptr<Parameters> params,
            const std::string outputFile,
            const int filename_count,
            const bool save_init,
            const bool onestep = false
            );
        void setStructure(std::shared_ptr<VPM::Structure> structure);

    private:

        bool new_stepsize(
                const double final_time,
                ParticleField & pf,
                double & delta_t
                );

        Point2d getForces(
                const ParticleField & pf,
                const ParticleField & pf_old,
                const double & delta_t
                );


        int m_dim;

        std::shared_ptr<Split_Diffusion> m_diffusion;
        std::shared_ptr<Split_Advection> m_advection;
        std::shared_ptr<Split_Source> m_source;

        bool m_structure_set;
        std::shared_ptr<VPM::Structure> m_structure;

        ParticleField m_pf;
        ParticleField m_pf_old;

        std::string m_outputFile;
        int m_filename_count;
        bool m_save_init;

    };

}
