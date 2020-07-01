#include "VPM2d.hpp"
#include "Point2d.hpp"
#include "InitialConditions.hpp"
#include "Parameters.hpp"
//#include "Structure_Flat.hpp"
//#include "Structure_Hilly.hpp"
//#include "Structure_Circle.hpp"
//#include "Structure_HalfCircle.hpp"
//#include "Structure_InverseCircle.hpp"
#include "Structure_Ellipse.hpp"
#include <string>

int main(int argc, char** argv)
{

    //std::string paramsFile;
    std::string inputFile;
    std::string outputFile;
    double final_time = 1;
    VPM::Point2d origo;
    double semimajoraxis;
    double semiminoraxis;
    bool random_velocity_dist = false;

    for(int i=1; i<argc; ) {
        int eat = 0;
        std::string arg( argv[i] );
        if( arg == "--if" && (i+1) < argc ) {
            inputFile = argv[i+1];
            eat = 2;
        }
        if( arg == "--rv" ) {
            random_velocity_dist = true;
            eat = 1;
        }
        //if( arg == "--pif" && (i+1) < argc ) {
        //    paramsFile = argv[i+1];
        //    eat = 2;
        //}
        if( arg == "--of" && (i+1) < argc ) {
            outputFile = argv[i+1];
            eat = 2;
        }
        if( arg == "--T" && (i+1) < argc ) {
            final_time = std::atof(argv[i+1]);
            eat = 2;
        }
        if( arg == "--origo" && (i+1) < argc ) {
            origo = VPM::Point2d(std::atof(argv[i+1]), std::atof(argv[i+2]));
            eat = 3;
        }
        if( arg == "--semimajoraxis" && (i+1) < argc ) {
            semimajoraxis = std::atof(argv[i+1]);
            eat = 2;
        }
        if( arg == "--semiminoraxis" && (i+1) < argc ) {
            semiminoraxis = std::atof(argv[i+1]);
            eat = 2;
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

    VPM::VPM2d* vpm = new VPM::VPM2d(argc, argv);

    std::shared_ptr<VPM::Structure_Ellipse> structure = std::make_shared<VPM::Structure_Ellipse>(origo, semimajoraxis, semiminoraxis);
    vpm->setStructure(structure);

    // ----------------- //
    // parameters //
    // ----------------- //
    //
    //std::shared_ptr<VPM::Parameters> params = std::make_shared<VPM::Parameters>();
    //params->readFromFile(paramsFile);

    // ----------------- //
    // position and vorticity //
    // ----------------- //
    //std::vector<VPM::Point2d> positions;
    //std::vector<double> omega;

    int fn_count;
    bool save_init;

    VPM::ParticleField pf;
    readParticlesFromFile(inputFile, pf, random_velocity_dist);//params->m_N, positions, omega);

    std::size_t found = inputFile.find_last_of("_");
    std::string num = inputFile.substr(found+2, inputFile.size());
    fn_count = std::atoi(num.c_str());
    save_init = false;

    // ----------------- //
    // run simulation //
    // ----------------- //
    vpm->run(pf, //positions, omega,
            final_time, //params,
            outputFile, fn_count, save_init, true
            );

    return 0;

}
