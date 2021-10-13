#include <optional>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <sstream>

#include "InitialConditions.hpp"
#include "Parameters.hpp"
#include "Particles2d.hpp"
#include "Point2d.hpp"
#include "Split_Advection.hpp"
#include "Structure.hpp"
#include "Structure_Ellipse.hpp"
#include "VPM2d.hpp"
#include <string>
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include <petsc.h>

std::string to_string(const VPM::Point2d &x) {
  return "Point2d(" + std::to_string(x.x) + ", " + std::to_string(x.y) + ")";
}

std::string to_string(const VPM::IPoint2d &x) {
  return "IPoint2d(" + std::to_string(x.x) + ", " + std::to_string(x.y) + ")";
}

PYBIND11_MODULE(pyVPM, m) {
  m.doc() = "Python Bindinds for the VPM module.";

  pybind11::class_<VPM::Point2d>(m, "Point2d")
      .def(pybind11::init<double, double>())
      .def_readwrite("x", &VPM::Point2d::x)
      .def_readwrite("y", &VPM::Point2d::y)
      .def("__repr__", [](const VPM::Point2d &x) { return to_string(x); })
      .def("tolist", [](const VPM::Point2d &x) {
        return std::vector<double>{{x.x, x.y}};
      });

  pybind11::class_<VPM::IPoint2d>(m, "IPoint2d")
      .def(pybind11::init<int, int>())
      .def_readwrite("x", &VPM::IPoint2d::x)
      .def_readwrite("y", &VPM::IPoint2d::y)
      .def("__repr__", [](const VPM::IPoint2d &x) { return to_string(x); })
      .def("tolist", [](const VPM::IPoint2d &x) {
        return std::vector<int>{{x.x, x.y}};
      });

  pybind11::class_<VPM::RemeshParams, std::shared_ptr<VPM::RemeshParams>>(
      m, "RemeshParams")
      .def(pybind11::init([](bool isOn, int steps, std::string method) {
             return std::make_shared<VPM::RemeshParams>(isOn, steps, method);
           }),
           pybind11::arg("isOn") = true, pybind11::arg("steps") = 1,
           pybind11::arg("method") = "M4d_fast");

  pybind11::class_<VPM::BCParams, std::shared_ptr<VPM::BCParams>>(m, "BCParams")
      .def(pybind11::init([](const std::string &xlstring, const double to_xl,
                             const std::string &xrstring, const double to_xr,
                             const std::string &ylstring, const double to_yl,
                             const std::string &yrstring, const double to_yr,
                             const std::string &zlstring, const double to_zl,
                             const std::string &zrstring, const double to_zr) {
             return std::make_shared<VPM::BCParams>(
                 xlstring, to_xl, xrstring, to_xr, ylstring, to_yl, yrstring,
                 to_yr, zlstring, to_zl, zrstring, to_zr);
           }),
           pybind11::arg("xlstring") = "none", pybind11::arg("to_xl") = 0.0,
           pybind11::arg("xrstring") = "none", pybind11::arg("to_xr") = 0.0,
           pybind11::arg("ylstring") = "none", pybind11::arg("to_yl") = 0.0,
           pybind11::arg("yrstring") = "none", pybind11::arg("to_yr") = 0.0,
           pybind11::arg("zlstring") = "none", pybind11::arg("to_zl") = 0.0,
           pybind11::arg("zrstring") = "none", pybind11::arg("to_zr") = 0.0);

  pybind11::class_<VPM::Parameters>(m, "Parameters")
      .def(pybind11::init(
          [](const VPM::Point2d &domain_ll, const VPM::Point2d &domain_ur,
             const unsigned int num_px, const unsigned int num_py,
             const double nu, const double population_threshold,
             const std::shared_ptr<VPM::RemeshParams> &remesh,
             const std::shared_ptr<VPM::BCParams> &bc,
             const unsigned int order_ODEsolver, const VPM::Point2d &Uinfty) {
            return VPM::Parameters(domain_ll, domain_ur, num_px, num_py, nu,
                                   population_threshold, remesh, bc,
                                   order_ODEsolver, Uinfty);
          }))
          .def("domain", [](const VPM::Parameters& params) {
              return std::vector<double>({params.m_domain_ll.x, params.m_domain_ur.x,
                params.m_domain_ll.x, params.m_domain_ur.y});
          })
          .def("nu", [](const VPM::Parameters& params) {
              return params.m_nu;
          })
          .def("nx", [](const VPM::Parameters& params) {
              return params.m_num_px; // FIXME: Is this the correct number?
          })
          .def("ny", [](const VPM::Parameters& params) {
              return params.m_num_py; // FIXME: Is this the correct number?
          });

  pybind11::class_<VPM::ParticleField, std::shared_ptr<VPM::ParticleField>>(m, "ParticleField")
      .def(pybind11::init([](){
        return std::make_shared<VPM::ParticleField>();
      }))
      .def_readwrite("positions", &VPM::ParticleField::positions)
      .def_readwrite("regular_positions",
                     &VPM::ParticleField::regular_positions)
      .def_readwrite("params", &VPM::ParticleField::params)
      .def_readwrite("omega", &VPM::ParticleField::omega)
      .def_readwrite("positions", &VPM::ParticleField::positions)
      .def_readwrite("velocity", &VPM::ParticleField::velocity)
      .def_readwrite("cartesianGrid", &VPM::ParticleField::cartesianGrid)
      .def_readwrite("velocity_correspondsTo_omega",
                     &VPM::ParticleField::velocity_correspondsTo_omega)
      .def_readwrite("Linfty_gradVelocity",
                     &VPM::ParticleField::Linfty_gradVelocity)
      .def_readwrite("time", &VPM::ParticleField::time);

  m.def("init_particles", &VPM::init, pybind11::arg("params"),
        pybind11::arg("example_num"), pybind11::arg("dist"),
        pybind11::arg("strength"), pybind11::arg("core_radius"),
        pybind11::arg("positions"), pybind11::arg("omega"));

  pybind11::class_<VPM::VPM2d>(m, "VPM2d")
      .def(pybind11::init([](const std::vector<std::string> &args) {
        // convert to char-vector
        std::vector<char *> temporary;
        for (const auto &arg : args) {
          char *data = const_cast<char *>(arg.c_str());
          temporary.push_back(data);
        }
        return VPM::VPM2d(temporary.size(), temporary.data());
      }))
      .def("run", &VPM::VPM2d::run, pybind11::arg("pf"),
           pybind11::arg("final_time"), pybind11::arg("outputFile"),
           pybind11::arg("fn_count"), pybind11::arg("save_init"),
           pybind11::arg("onestep"))
      .def("run_without_writer", &VPM::VPM2d::run_without_writer, pybind11::arg("pf"),
           pybind11::arg("final_time"),
           pybind11::arg("onestep"))
      .def("setStructure", &VPM::VPM2d::setStructure,
           pybind11::arg("setStructure"));

  pybind11::class_<VPM::Split_Advection>(m, "Split_Advection")
      .def("calculateVelocity", &VPM::Split_Advection::calculateVelocity,
           pybind11::arg("pf"));

  pybind11::class_<VPM::Structure, std::shared_ptr<VPM::Structure>>(
      m, "Structure");
  pybind11::class_<VPM::Structure_Ellipse,
                   std::shared_ptr<VPM::Structure_Ellipse>, VPM::Structure>(
      m, "Structure_Ellipse")
      .def(pybind11::init([](VPM::Point2d origo, double semi_major_axis,
                             double semi_minor_axis) {
             return std::make_shared<VPM::Structure_Ellipse>(
                 origo, semi_major_axis, semi_minor_axis);
           }),
           pybind11::arg("origo"), pybind11::arg("semi_major_axis"),
           pybind11::arg("semi_minor_axis"));

  m.def("initialize_module", [](const std::vector<std::string> &args) {
    // convert to char-vector
    std::vector<char *> temporary;
    for (const auto &arg : args) {
      char *data = const_cast<char *>(arg.c_str());
      temporary.push_back(data);
    }
    int new_size = temporary.size();
    char **temporary_data = temporary.data();
    MPI_Init(&new_size, &temporary_data);
  });

  m.def("readParticlesFromFile", VPM::readParticlesFromFile,
        pybind11::arg("filename"), pybind11::arg("pf"),
        pybind11::arg("random_velocity_dist"));
}
