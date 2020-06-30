#include "VPM3d.hpp"
#include "Point3d.hpp"
#include "InitialConditions.hpp"
#include "Parameters.hpp"
#include "Structure_Flat.hpp"
#include "Structure_Hilly.hpp"
#include "Structure_Sphere.hpp"
#include "Structure_Bar.hpp"
#include <string>
#include <memory>

int m_numberofParticlesx = 50;
int m_numberofParticlesy = 50;
int m_numberofParticlesz = 50;

int main(int argc, char** argv)
{

    std::string inputFile;
    std::string outputFile;
    double time = 0.1;
    double dt = 0.01;
    double nu = 0.;

    double ax = 1.0;
    double ay = 1.0;
    double az = 1.0;
    VPM::Point3d domain_ll = VPM::Point3d(-ax,-ay,-az);
    VPM::Point3d domain_ur = VPM::Point3d(+ax,+ay,+az);

    VPM::Point3d Uinfty = VPM::Point3d(0.);

    bool remesh_isOn = true;
    int remesh_steps = 1;
    std::string remesh_method = "order6_fast";

    int order = 1;

    int example_num = 0;
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
    std::string bc_zl = "none";
    double bc_to_zl = double(0.);
    std::string bc_zr = "none";
    double bc_to_zr = double(0.);

    for(int i=1; i<argc; ) {
        int eat = 0;
        std::string arg( argv[i] );
        if( arg == "--if" && (i+1) < argc )
        {
            inputFile = argv[i+1];
            eat = 2;
        }
        if( arg == "--of" && (i+1) < argc ) {
            outputFile = argv[i+1];
            eat = 2;
        }
        if( arg == "--remesh_method" && (i+1) < argc ) {
            remesh_method = argv[i+1];
            eat = 2;
        }
        if( arg == "--population_threshold" && (i+1) < argc ) {
            population_threshold = std::atof(argv[i+1]);
            eat = 2;
        }
        if( arg == "--nu" && (i+1) < argc ) {
            nu = std::atof(argv[i+1]);
            eat = 2;
        }
        if( arg == "--T" && (i+1) < argc ) {
            time = std::atof(argv[i+1]);
            eat = 2;
        }
        //if( arg == "--dt" && (i+1) < argc ) {
        //    dt = std::atof(argv[i+1]);
        //    eat = 2;
        //}
        if( arg == "--domain" && (i+6) < argc ) {
            domain_ll = VPM::Point3d(std::atof(argv[i+1]),std::atof(argv[i+2]),std::atof(argv[i+3]));
            domain_ur = VPM::Point3d(std::atof(argv[i+4]),std::atof(argv[i+5]),std::atof(argv[i+6]));
            eat = 7;
        }
        if( arg == "--Uinfty" && (i+3) < argc ) {
            Uinfty = VPM::Point3d(std::atof(argv[i+1]),std::atof(argv[i+2]),std::atof(argv[i+3]));
            eat = 4;
        }
        if( arg == "--nx" && (i+1) < argc ) {
            m_numberofParticlesx = std::atoi(argv[i+1]);
            eat = 2;
        }
        if( arg == "--ny" && (i+1) < argc ) {
            m_numberofParticlesy = std::atoi(argv[i+1]);
            eat = 2;
        }
        if( arg == "--nz" && (i+1) < argc ) {
            m_numberofParticlesz = std::atoi(argv[i+1]);
            eat = 2;
        }
        if( arg == "--order" && (i+1) < argc ) {
            order = std::atoi(argv[i+1]);
            eat = 2;
        }
        if( arg == "--example_num" && (i+1) < argc ) {
            example_num = std::atoi(argv[i+1]);
            eat = 2;
        }
        if( arg == "--structure_num" && (i+1) < argc ) {
            structure_num = std::atoi(argv[i+1]);
            eat = 2;
        }
        if( arg == "--example_dist" && (i+1) < argc ) {
            example_dist = std::atof(argv[i+1]);
            eat = 2;
        }
        if( arg == "--example_strength" && (i+1) < argc ) {
            example_strength = std::atof(argv[i+1]);
            eat = 2;
        }
        if( arg == "--example_core_radius" && (i+1) < argc ) {
            example_core_radius = std::atof(argv[i+1]);
            eat = 2;
        }
        if( arg == "--bc_xl" && (i+2) < argc ) {
            bc_xl = argv[i+1];
            bc_to_xl = std::atof(argv[i+2]);
            eat = 3;
        }
        if( arg == "--bc_xr" && (i+2) < argc ) {
            bc_xr = argv[i+1];
            bc_to_xr = std::atof(argv[i+2]);
            eat = 3;
        }
        if( arg == "--bc_yl" && (i+2) < argc ) {
            bc_yl = argv[i+1];
            bc_to_yl = std::atof(argv[i+2]);
            eat = 3;
        }
        if( arg == "--bc_yr" && (i+2) < argc ) {
            bc_yr = argv[i+1];
            bc_to_yr = std::atof(argv[i+2]);
            eat = 3;
        }
        if( arg == "--bc_zl" && (i+2) < argc ) {
            bc_zl = argv[i+1];
            bc_to_zl = std::atof(argv[i+2]);
            eat = 3;
        }
        if( arg == "--bc_zr" && (i+2) < argc ) {
            bc_zr = argv[i+1];
            bc_to_zr = std::atof(argv[i+2]);
            eat = 3;
        }
        if( eat > 0 ) {
            argc = argc - eat;
            for(int k=i; k<argc; k++ ) {
                argv[k] = argv[k+eat];
            }
        }
        else {
            i++;
        }
    }

    VPM::VPM3d* vpm = new VPM::VPM3d(argc, argv);

    switch (structure_num)
    {
        case (0):
            {
                std::shared_ptr<VPM::Structure_Flat> structure = std::make_shared<VPM::Structure_Flat>();
                vpm->setStructure(structure);
            }
            break;
        case (1):
            {
                std::shared_ptr<VPM::Structure_Hilly> structure = std::make_shared<VPM::Structure_Hilly>();
                vpm->setStructure(structure);
            }
            break;
        case (2):
            {
                std::shared_ptr<VPM::Structure_Sphere> structure = std::make_shared<VPM::Structure_Sphere>();
                vpm->setStructure(structure);
            }
            break;
        case (3):
            {
                std::shared_ptr<VPM::Structure_Bar> structure = std::make_shared<VPM::Structure_Bar>();
                vpm->setStructure(structure);
            }
            break;
    }

    std::shared_ptr<VPM::RemeshParams> remeshParams = std::make_shared<VPM::RemeshParams>(
            remesh_isOn,
            remesh_steps,
            remesh_method
            );

    std::shared_ptr<VPM::BCParams> bcParams = std::make_shared<VPM::BCParams>(
            bc_xl,
            bc_to_xl,
            bc_xr,
            bc_to_xr,
            bc_yl,
            bc_to_yl,
            bc_yr,
            bc_to_yr,
            bc_zl,
            bc_to_zl,
            bc_zr,
            bc_to_zr
            );

    std::shared_ptr<VPM::Parameters> params;
    std::vector<VPM::Point3d> positions;
    std::vector<VPM::Point3d> omega;

        params = std::make_shared<VPM::Parameters>(
                domain_ll,
                domain_ur,
                m_numberofParticlesx,
                m_numberofParticlesy,
                m_numberofParticlesz,
                nu,
                0,
                population_threshold,
                remeshParams,
                bcParams,
                order,
                Uinfty
                );

    int fn_count;
    bool save_init;
    if (inputFile.empty())
    {
        init(params, example_num, example_dist, example_strength, example_core_radius, positions, omega);
        fn_count = 0;
        save_init = true;

    }
    else
    {
        params->readFromFile(inputFile+"_conf.dat");
        readParticlesFromFile(inputFile, params, positions, omega);

        std::size_t found = inputFile.find_last_of("_");
        std::string num = inputFile.substr(found+2, inputFile.size());
        fn_count = std::atoi(num.c_str());
        save_init = false;
    }

    vpm->run(positions, omega, time, dt, params,
            outputFile, fn_count, save_init
            );

    return 0;

}
