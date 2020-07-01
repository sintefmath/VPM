#include "InitialConditions.hpp"
#include "FlattenArrayIndex.hpp"
#include <iostream>
#include <fstream>

namespace VPM {

    /*
     * omega(r,t) = Gamma_0/(pi*r_c(t)^2) * exp(-r^2/r_c(t)^2)
     * v_theta(r,t) = Gamma_0/(2*pi*r) * (1 - exp(-r^2/r_c(t)^2))
     */
    double omega_LambOseen(const Point2d x, const double strength, const double core_radius)
    {
        double r2 = dot(x,x);
        //double t = 1.;
        //double nu = core_radius;
        //double rct_squared = 4.*nu*t;
        double rct_squared = 0.15*0.15;
        return strength/(M_PI*rct_squared) * exp(-r2/rct_squared);
    }


    double omega_CoRotatingVortexPair(const Point2d x, const double dist, const double strength, const double core_radius)
    {
        Point2d c_a = Point2d(-dist/2.,0.);
        Point2d c_b = Point2d(+dist/2.,0.);
        return omega_LambOseen(x-c_a, strength, core_radius) + omega_LambOseen(x-c_b, strength, core_radius);
    }
    //double omega_CoRotatingVortexPairSharp(const Point2d x, const double dist)
    //{
    //    Point2d c_a = Point2d(-dist/2.,0.);
    //    Point2d c_b = Point2d(+dist/2.,0.);
    //    double ra = dot(x-c_a,x-c_a);
    //    double rb = dot(x-c_b,x-c_b);
    //    double val = 0;
    //    if (ra<dist/12. || rb<dist/12.)
    //    {
    //        val = 1.0;
    //    }
    //    return val;
    //}
    double omega_CounterRotatingVortexPair(const Point2d x, const double dist, const double strength, const double core_radius)
    {
        //Point2d c_a = Point2d(-dist/2.,.0);
        //Point2d c_b = Point2d(+dist/2.,.0);
        Point2d c_a = Point2d(-dist/2.,+.5);
        Point2d c_b = Point2d(+dist/2.,+.5);
        return - omega_LambOseen(x-c_a, strength, core_radius) + omega_LambOseen(x-c_b, strength, core_radius);
    }

    void init(
            const VPM::Parameters& params,
            const int example_num,
            const double dist,
            const double strength,
            const double core_radius,
            std::vector<Point2d> & positions,
            std::vector<double> & omega
            )
    {

        positions.clear();
        omega.clear();

        positions.resize(params.m_N);
        omega.resize(params.m_N);

        int n = -1;
        for (unsigned int npy=0; npy<params.m_num_py; npy++)
        {
            for (unsigned int npx=0; npx<params.m_num_px; npx++)
            {
                Point2d pos = params.m_domain_ll + Point2d(
                        double(2*npx+1)/double(2*params.m_num_px)*(params.m_domain_ur.x-params.m_domain_ll.x),
                        double(2*npy+1)/double(2*params.m_num_py)*(params.m_domain_ur.y-params.m_domain_ll.y)
                        );
                double w;
                switch (example_num)
                {
                    case -1:
                        w = double(0.);
                        break;
                    case 0:
                        w = omega_LambOseen(pos, strength, core_radius);
                        break;
                    case 1:
                        w = omega_CoRotatingVortexPair(pos, dist, strength, core_radius);
                        break;
                    case 2:
                        w = omega_CounterRotatingVortexPair(pos, dist, strength, core_radius);
                        break;
                    //case 3:
                    //    w = omega_CoRotatingVortexPairSharp(pos, dist);
                    //    break;
                }
                positions[flat(npx,npy,int(params.m_num_px))] = pos;
                omega[flat(npx,npy,int(params.m_num_px))] = w;
            }
        }

    }

    //bool init(const std::string & filename)
    //{
    //
    //
    //    bool success = readParticlesFromFile(filename);
    //    if (success)
    //    {
    //        updateParams();
    //        return true;
    //    }
    //    else
    //    {
    //        return false;
    //    }
    //}

    void writeParticlesToFile(
            const std::string & filename,
            const int timestep,
            const ParticleField & pf
            )
    {
        if (!pf.velocity_correspondsTo_omega)
        {
            std::cerr<<"---> Warning: Velocity field does not correspond to vorticity"<<std::endl;
        }

        std::string t = std::to_string(timestep);

        pf.params.writeToFile(filename+"_t"+t+"_conf.dat");

        std::ofstream myfile_posx;
        std::ofstream myfile_posy;
        std::ofstream myfile_regposx;
        std::ofstream myfile_regposy;
        std::ofstream myfile_omega;
        std::ofstream myfile_velx;
        std::ofstream myfile_vely;
        std::ofstream myfile_cartesianGrid;
        std::ofstream myfile_velocity_correspondsTo_omega;
        std::ofstream myfile_Linfty_gradVelocity;
        std::ofstream myfile_pf_time;

        myfile_posx.open (filename+"_t"+t+"_posx.dat", std::ofstream::binary);
        myfile_posy.open (filename+"_t"+t+"_posy.dat", std::ofstream::binary);
        myfile_regposx.open (filename+"_t"+t+"_regposx.dat", std::ofstream::binary);
        myfile_regposy.open (filename+"_t"+t+"_regposy.dat", std::ofstream::binary);
        myfile_omega.open (filename+"_t"+t+"_omega.dat", std::ofstream::binary);
        myfile_velx.open (filename+"_t"+t+"_Ux.dat", std::ofstream::binary);
        myfile_vely.open (filename+"_t"+t+"_Uy.dat", std::ofstream::binary);
        myfile_cartesianGrid.open (filename+"_t"+t+"_cartesianGrid.dat", std::ofstream::binary);
        myfile_velocity_correspondsTo_omega.open (filename+"_t"+t+"_velocity_correspondsTo_omega.dat", std::ofstream::binary);
        myfile_Linfty_gradVelocity.open (filename+"_t"+t+"_Linfty_gradVelocity.dat", std::ofstream::binary);
        myfile_pf_time.open (filename+"_t"+t+"_time.dat", std::ofstream::binary);

        for (unsigned int i=0; i<pf.positions.size(); i++)
        {
            myfile_posx.write(reinterpret_cast<const char*>(&pf.positions[i].x), sizeof(double));
            myfile_posy.write(reinterpret_cast<const char*>(&pf.positions[i].y), sizeof(double));
            myfile_omega.write(reinterpret_cast<const char*>(&pf.omega[i]), sizeof(double));
            myfile_velx.write(reinterpret_cast<const char*>(&pf.velocity[i].x), sizeof(double));
            myfile_vely.write(reinterpret_cast<const char*>(&pf.velocity[i].y), sizeof(double));
        }

        for (unsigned int i=0; i<pf.regular_positions.size(); i++)
        {
            myfile_regposx.write(reinterpret_cast<const char*>(&pf.regular_positions[i].x), sizeof(double));
            myfile_regposy.write(reinterpret_cast<const char*>(&pf.regular_positions[i].y), sizeof(double));
        }

        myfile_cartesianGrid.write(reinterpret_cast<const char*>(&pf.cartesianGrid), sizeof(bool));
        myfile_velocity_correspondsTo_omega.write(reinterpret_cast<const char*>(&pf.velocity_correspondsTo_omega), sizeof(bool));
        myfile_Linfty_gradVelocity.write(reinterpret_cast<const char*>(&pf.Linfty_gradVelocity), sizeof(double));
        myfile_pf_time.write(reinterpret_cast<const char*>(&pf.time), sizeof(double));

        myfile_posx.close();
        myfile_posy.close();
        myfile_regposx.close();
        myfile_regposy.close();
        myfile_omega.close();
        myfile_velx.close();
        myfile_vely.close();
        myfile_cartesianGrid.close();
        myfile_velocity_correspondsTo_omega.close();
        myfile_Linfty_gradVelocity.close();
        myfile_pf_time.close();

    }

    bool readParticlesFromFile(
            const std::string & filename,
            ParticleField & pf,
            bool random_velocity_dist
            )
    {

        pf.params.readFromFile(filename+"_conf.dat");

        std::ifstream myfile_posx(filename+"_posx.dat", std::ifstream::binary);
        if(!myfile_posx.is_open())
        {
            std::cerr<<"\nCannot open file: "<< filename<<"_posx.dat";
            return false;
        }
        std::ifstream myfile_posy(filename+"_posy.dat", std::ifstream::binary);
        if(!myfile_posy.is_open())
        {
            std::cerr<<"\nCannot open file: "<< filename<<"_posy.dat";
            return false;
        }
        std::ifstream myfile_regposx(filename+"_regposx.dat", std::ifstream::binary);
        if(!myfile_regposx.is_open())
        {
            std::cerr<<"\nCannot open file: "<< filename<<"_regposx.dat";
            return false;
        }
        std::ifstream myfile_regposy(filename+"_regposy.dat", std::ifstream::binary);
        if(!myfile_regposy.is_open())
        {
            std::cerr<<"\nCannot open file: "<< filename<<"_regposy.dat";
            return false;
        }
        std::ifstream myfile_omega(filename+"_omega.dat", std::ifstream::binary);
        if(!myfile_omega.is_open())
        {
            std::cerr<<"\nCannot open file: "<< filename<<"_omega.dat";
            return false;
        }
        std::ifstream myfile_velx(filename+"_Ux.dat", std::ifstream::binary);
        if(!myfile_velx.is_open())
        {
            std::cerr<<"\nCannot open file: "<< filename<<"_Ux.dat";
            return false;
        }
        std::ifstream myfile_vely(filename+"_Uy.dat", std::ifstream::binary);
        if(!myfile_vely.is_open())
        {
            std::cerr<<"\nCannot open file: "<< filename<<"_Uy.dat";
            return false;
        }
        std::ifstream myfile_cartesianGrid(filename+"_cartesianGrid.dat", std::ifstream::binary);
        if(!myfile_cartesianGrid.is_open())
        {
            std::cerr<<"\nCannot open file: "<< filename<<"_cartesianGrid.dat";
            return false;
        }
        std::ifstream myfile_velocity_correspondsTo_omega(filename+"_velocity_correspondsTo_omega.dat", std::ifstream::binary);
        if(!myfile_velocity_correspondsTo_omega.is_open())
        {
            std::cerr<<"\nCannot open file: "<< filename<<"_velocity_correspondsTo_omega.dat";
            return false;
        }
        std::ifstream myfile_Linfty_gradVelocity(filename+"_Linfty_gradVelocity.dat", std::ifstream::binary);
        if(!myfile_Linfty_gradVelocity.is_open())
        {
            std::cerr<<"\nCannot open file: "<< filename<<"_Linfty_gradVelocity.dat";
            return false;
        }
        std::ifstream myfile_pf_time(filename+"_time.dat", std::ifstream::binary);
        if(!myfile_pf_time.is_open())
        {
            std::cerr<<"\nCannot open file: "<< filename<<"_pf_time.dat";
            return false;
        }

        pf.positions.clear();
        pf.regular_positions.clear();
        pf.omega.clear();
        pf.velocity.clear();

        double ra=0.;
        for (unsigned int i=0; i<pf.params.m_N; i++)
        {
            double posx, posy;
            double tmp_omega;
            double velx, vely;

            myfile_posx.read((char *) &posx, sizeof(double));
            myfile_posy.read((char *) &posy, sizeof(double));
            myfile_omega.read((char *) &tmp_omega, sizeof(double));
            myfile_velx.read((char *) &velx, sizeof(double));
            myfile_vely.read((char *) &vely, sizeof(double));

            if (random_velocity_dist)
            {
                ra = (2*((double) rand()/RAND_MAX)-1.)/1.;
            }

            pf.positions.push_back(Point2d(posx, posy));
            pf.omega.push_back(tmp_omega+ra);
            pf.velocity.push_back(Point2d(velx, vely));
        }

        for (unsigned int i=0; i<pf.params.m_num_px*pf.params.m_num_py; i++)
        {
            double reg_posx, reg_posy;

            myfile_regposx.read((char *) &reg_posx, sizeof(double));
            myfile_regposy.read((char *) &reg_posy, sizeof(double));
            pf.regular_positions.push_back(Point2d(reg_posx, reg_posy));
        }

        myfile_cartesianGrid.read((char *) &pf.cartesianGrid, sizeof(bool));
        if (random_velocity_dist)
        {
            pf.velocity_correspondsTo_omega = false;
        }
        else
        {
            myfile_velocity_correspondsTo_omega.read((char *) &pf.velocity_correspondsTo_omega, sizeof(bool));
        }
        myfile_Linfty_gradVelocity.read((char *) &pf.Linfty_gradVelocity, sizeof(double));
        myfile_pf_time.read((char *) &pf.time, sizeof(double));

        myfile_posx.close();
        myfile_posy.close();
        myfile_omega.close();
        myfile_regposx.close();
        myfile_regposy.close();
        myfile_cartesianGrid.close();
        myfile_velocity_correspondsTo_omega.close();
        myfile_Linfty_gradVelocity.close();
        myfile_pf_time.close();
        return true;
    }

}
