#include "InitialConditions.hpp"
#include "Parameters.hpp"
#include "Point2d.hpp"
#include "Structure_Circle.hpp"
#include "Structure_Ellipse.hpp"
#include "Structure_Flat.hpp"
#include "Structure_HalfCircle.hpp"
#include "Structure_Hilly.hpp"
#include "Structure_InverseCircle.hpp"
#include "VPM2d.hpp"
#define OMPI_SKIP_MPICXX 1

#include <mpi.h>
#include <string>

int m_numberofParticlesx = 50;
int m_numberofParticlesy = 50;

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  std::string inputFile;
  std::string outputFile;
  double time = 0.1;
  double dt = 0.01;
  double nu = 0.;

  double ax = 1.0;
  double ay = 1.0;
  VPM::Point2d domain_ll = VPM::Point2d(-ax, -ay);
  VPM::Point2d domain_ur = VPM::Point2d(+ax, +ay);

  VPM::Point2d Uinfty = VPM::Point2d(0.);

  bool remesh_isOn = true;
  int remesh_steps = 1;
  std::string remesh_method = "M6d_fast";

  int order = 1;

  int example_num = -1;
  int structure_num = -1;
  double example_dist = 0.5;
  double example_strength = 1.;
  double example_core_radius = 0.15;

  double population_threshold = -1;

  std::string bc_xl = "none";
  double bc_to_xl = double(0.);
  std::string bc_xr = "none";
  double bc_to_xr = double(0.);
  std::string bc_yl = "none";
  double bc_to_yl = double(0.);
  std::string bc_yr = "none";
  double bc_to_yr = double(0.);

  for (int i = 1; i < argc;) {
    int eat = 0;
    std::string arg(argv[i]);
    if (arg == "--if" && (i + 1) < argc) {
      inputFile = argv[i + 1];
      eat = 2;
    }
    if (arg == "--of" && (i + 1) < argc) {
      outputFile = argv[i + 1];
      eat = 2;
    }
    if (arg == "--remesh_method" && (i + 1) < argc) {
      remesh_method = argv[i + 1];
      eat = 2;
    }
    if (arg == "--population_threshold" && (i + 1) < argc) {
      population_threshold = std::atof(argv[i + 1]);
      eat = 2;
    }
    if (arg == "--nu" && (i + 1) < argc) {
      nu = std::atof(argv[i + 1]);
      eat = 2;
    }
    if (arg == "--T" && (i + 1) < argc) {
      time = std::atof(argv[i + 1]);
      eat = 2;
    }
    if (arg == "--dt" && (i + 1) < argc) {
      dt = std::atof(argv[i + 1]);
      eat = 2;
    }
    if (arg == "--domain" && (i + 4) < argc) {
      domain_ll = VPM::Point2d(std::atof(argv[i + 1]), std::atof(argv[i + 2]));
      domain_ur = VPM::Point2d(std::atof(argv[i + 3]), std::atof(argv[i + 4]));
      eat = 5;
    }
    if (arg == "--Uinfty" && (i + 2) < argc) {
      Uinfty = VPM::Point2d(std::atof(argv[i + 1]), std::atof(argv[i + 2]));
      eat = 3;
    }
    if (arg == "--nx" && (i + 1) < argc) {
      m_numberofParticlesx = std::atoi(argv[i + 1]);
      eat = 2;
    }
    if (arg == "--ny" && (i + 1) < argc) {
      m_numberofParticlesy = std::atoi(argv[i + 1]);
      eat = 2;
    }
    if (arg == "--order" && (i + 1) < argc) {
      order = std::atoi(argv[i + 1]);
      eat = 2;
    }
    if (arg == "--example_num" && (i + 1) < argc) {
      example_num = std::atoi(argv[i + 1]);
      eat = 2;
    }
    if (arg == "--structure_num" && (i + 1) < argc) {
      structure_num = std::atoi(argv[i + 1]);
      eat = 2;
    }
    if (arg == "--example_dist" && (i + 1) < argc) {
      example_dist = std::atof(argv[i + 1]);
      eat = 2;
    }
    if (arg == "--example_strength" && (i + 1) < argc) {
      example_strength = std::atof(argv[i + 1]);
      eat = 2;
    }
    if (arg == "--example_core_radius" && (i + 1) < argc) {
      example_core_radius = std::atof(argv[i + 1]);
      eat = 2;
    }
    if (arg == "--bc_xl" && (i + 2) < argc) {
      bc_xl = argv[i + 1];
      bc_to_xl = std::atof(argv[i + 2]);
      eat = 3;
    }
    if (arg == "--bc_xr" && (i + 2) < argc) {
      bc_xr = argv[i + 1];
      bc_to_xr = std::atof(argv[i + 2]);
      eat = 3;
    }
    if (arg == "--bc_yl" && (i + 2) < argc) {
      bc_yl = argv[i + 1];
      bc_to_yl = std::atof(argv[i + 2]);
      eat = 3;
    }
    if (arg == "--bc_yr" && (i + 2) < argc) {
      bc_yr = argv[i + 1];
      bc_to_yr = std::atof(argv[i + 2]);
      eat = 3;
    }
    if (eat > 0) {
      argc = argc - eat;
      for (int k = i; k < argc; k++) {
        argv[k] = argv[k + eat];
      }
    } else {
      i++;
    }
  }

  VPM::VPM2d *vpm = new VPM::VPM2d(argc, argv);

  switch (structure_num) {
  case (0): {
    std::shared_ptr<VPM::Structure_Flat> structure =
        std::make_shared<VPM::Structure_Flat>();
    vpm->setStructure(structure);
  } break;
  case (1): {
    std::shared_ptr<VPM::Structure_Hilly> structure =
        std::make_shared<VPM::Structure_Hilly>();
    vpm->setStructure(structure);
  } break;
  case (2): {
    std::shared_ptr<VPM::Structure_Circle> structure =
        std::make_shared<VPM::Structure_Circle>();
    vpm->setStructure(structure);
  } break;
  case (3): {
    std::shared_ptr<VPM::Structure_HalfCircle> structure =
        std::make_shared<VPM::Structure_HalfCircle>();
    vpm->setStructure(structure);
  } break;
  case (4): {
    std::shared_ptr<VPM::Structure_InverseCircle> structure =
        std::make_shared<VPM::Structure_InverseCircle>();
    vpm->setStructure(structure);
  } break;
  case (5): {
    std::shared_ptr<VPM::Structure_Ellipse> structure =
        std::make_shared<VPM::Structure_Ellipse>();
    vpm->setStructure(structure);
  } break;
  }

  std::shared_ptr<VPM::RemeshParams> remeshParams =
      std::make_shared<VPM::RemeshParams>(remesh_isOn, remesh_steps,
                                          remesh_method);

  std::shared_ptr<VPM::BCParams> bcParams = std::make_shared<VPM::BCParams>(
      bc_xl, bc_to_xl, bc_xr, bc_to_xr, bc_yl, bc_to_yl, bc_yr, bc_to_yr);

  std::shared_ptr<VPM::Parameters> params = std::make_shared<VPM::Parameters>(
      domain_ll, domain_ur, m_numberofParticlesx, m_numberofParticlesy, nu,
      // 0, // FIXME: This seems not to be used?
      population_threshold, remeshParams, bcParams, order, Uinfty);

  std::vector<VPM::Point2d> positions;
  std::vector<double> omega;
  VPM::ParticleField pf;
  int fn_count;
  bool save_init;
  if (inputFile.empty()) {
    init(*params, example_num, example_dist, example_strength,
         example_core_radius, positions, omega);
    pf.positions = positions;
    fn_count = 0;
    save_init = true;

  } else {

    readParticlesFromFile(inputFile, pf);

    std::size_t found = inputFile.find_last_of("_");
    std::string num = inputFile.substr(found + 2, inputFile.size());
    fn_count = std::atoi(num.c_str());
    save_init = false;
  }

  vpm->run(pf,
           // omega, // FIXME: This seems unused
           time,
           // dt, // FIXME: This seems unused
           // params, // FIXME: This seems unused
           outputFile, fn_count, save_init);

  return 0;
}
