#pragma once
#include "Particles3d.hpp"
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
    class VPM3d
    {
    public:
        VPM3d(int argc, char** argv);
        ~VPM3d();

        void run(
            const std::vector<Point3d> & positions,
            const std::vector<Point3d> & omega,
            const double T,
            const double dt,
            const std::shared_ptr<Parameters> params,
            const std::string outputFile,
            const int filename_count,
            const bool save_init
            );

        void setStructure(std::shared_ptr<VPM::Structure> structure);

    private:

        bool new_stepsize(
                const double T,
                ParticleField & pf,
                double & delta_t
                );

        int m_dim;

        std::shared_ptr<Split_Diffusion> m_diffusion;
        std::shared_ptr<Split_Advection> m_advection;
        std::shared_ptr<Split_Source> m_source;

        ParticleField m_pf;

        std::string m_outputFile;
        int m_filename_count;
        bool m_save_init;

    };

}
