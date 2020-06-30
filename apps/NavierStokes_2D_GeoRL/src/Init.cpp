#include "VPM2d.hpp"
#include "Particles2d.hpp"
#include "Point2d.hpp"
#include "Split_Advection.hpp"
#include "Parameters.hpp"
#include <string>

int main(int argc, char** argv)
{

    std::string outputFile = "init.dat";

    double dt = 0.01;
    double nu = 0.;
    int m_numberofParticlesx = 50;
    int m_numberofParticlesy = 50;


    double ax = 1.0;
    double ay = 1.0;
    VPM::Point2d domain_ll = VPM::Point2d(-ax,-ay);
    VPM::Point2d domain_ur = VPM::Point2d(+ax,+ay);

    VPM::Point2d Uinfty = VPM::Point2d(0.);

    bool remesh_isOn = true;
    int remesh_steps = 1;
    std::string remesh_method = "M6d_fast";

    int order = 1;

    double population_threshold = -1;

    std::string bc_xl = "none";
    double bc_to_xl = double(0.);
    std::string bc_xr = "none";
    double bc_to_xr = double(0.);
    std::string bc_yl = "none";
    double bc_to_yl = double(0.);
    std::string bc_yr = "none";
    double bc_to_yr = double(0.);

    for(int i=1; i<argc; ) {
        int eat = 0;
        std::string arg( argv[i] );
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
        if( arg == "--domain" && (i+4) < argc ) {
            domain_ll = VPM::Point2d(std::atof(argv[i+1]),std::atof(argv[i+2]));
            domain_ur = VPM::Point2d(std::atof(argv[i+3]),std::atof(argv[i+4]));
            eat = 5;
        }
        if( arg == "--Uinfty" && (i+2) < argc ) {
            Uinfty = VPM::Point2d(std::atof(argv[i+1]),std::atof(argv[i+2]));
            eat = 3;
        }
        if( arg == "--nx" && (i+1) < argc ) {
            m_numberofParticlesx = std::atoi(argv[i+1]);
            eat = 2;
        }
        if( arg == "--ny" && (i+1) < argc ) {
            m_numberofParticlesy = std::atoi(argv[i+1]);
            eat = 2;
        }
        if( arg == "--order" && (i+1) < argc ) {
            order = std::atoi(argv[i+1]);
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
            bc_to_yr
            );

    VPM::Parameters params = VPM::Parameters(
            domain_ll,
            domain_ur,
            m_numberofParticlesx,
            m_numberofParticlesy,
            nu,
            population_threshold,
            remeshParams,
            bcParams,
            order,
            Uinfty
            );

    //params->writeToFile(paramsFile);
    //params.print();

    std::vector<VPM::Point2d> positions;
    std::vector<double> omega;
    init(params, -1, -1, -1, -1, positions, omega);

    VPM::ParticleField pf;
    pf.omega = omega;
    pf.positions = positions;
    pf.regular_positions = positions;
    pf.params = params;
    pf.cartesianGrid = true;
    pf.velocity_correspondsTo_omega = false;

    VPM::Split_Advection advection;
    advection.calculateVelocity(pf);

    // it is important to set this, because a structure can/will be set
    pf.velocity_correspondsTo_omega = false;
    std::cerr<<"velocity_correspondsTo_omega="<<pf.velocity_correspondsTo_omega;

    writeParticlesToFile(outputFile, 0, pf);

    return 0;

}
