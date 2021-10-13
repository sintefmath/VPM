#pragma once
#include "InitialConditions.hpp"
#include "Parameters.hpp"
#include "Particles2d.hpp"
#include "Structure.hpp"

#include <memory>
#include <vector>
#include <functional>

namespace VPM {
class Split_Advection;
class Split_Diffusion;
class Split_Source;
} // namespace VPM

namespace VPM {
class VPM2d {
public:
  VPM2d(int argc, char **argv);
  ~VPM2d();

  void run(ParticleField &pf,
           // const std::vector<Point2d> & positions,
           // const std::vector<double> & omega,
           const double final_time,
           // const std::shared_ptr<Parameters> params,
           const std::string outputFile, const int filename_count,
           const bool save_init, const bool onestep = false);

  void run_with_writer(ParticleField &pf, const double final_time,
                  const bool onestep,
                  std::function<void(int, const ParticleField &)> writer);

  void run_without_writer(ParticleField &pf, const double final_time,
                  const bool onestep);                  
  void setStructure(std::shared_ptr<VPM::Structure> structure);

private:
  bool new_stepsize(const double final_time, ParticleField &pf,
                    double &delta_t);

  int m_dim;

  std::shared_ptr<Split_Diffusion> m_diffusion;
  std::shared_ptr<Split_Advection> m_advection;
  std::shared_ptr<Split_Source> m_source;

  ParticleField m_pf;

  std::string m_outputFile;
  int m_filename_count = 0;
  bool m_save_init = false;
};

} // namespace VPM
