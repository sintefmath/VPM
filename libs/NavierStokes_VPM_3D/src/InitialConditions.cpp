#include "InitialConditions.hpp"
#include <iostream>
#include <fstream>

namespace VPM {

    //Point3d cartesian2polar(const Point3d x)
    //{
    //    double r = L2Norm(x);
    //    if (r<1e-16)
    //    {
    //        return Point3d(0,0,0);
    //    }
    //    else
    //    {
    //        // x ~ r
    //        // y ~ theta = azimuth
    //        // z ~ phi = polar
    //        return Point3d(
    //                r,
    //                atan2(x.y, x.x),
    //                acos(x.z/r)
    //                );
    //    }
    //}
    //Point3d polar2cartesian(const Point3d x)
    //{
    //    // x ~ r
    //    // y ~ theta = azimuth
    //    // z ~ phi = polar
    //    return Point3d(
    //            x.x*cos(x.y)*sin(x.z),
    //            x.x*sin(x.y)*sin(x.z), 
    //            x.x*         cos(x.z)
    //            );
    //}

    Point3d omega_VortexRing_z(const Point3d x, const double strength, const double R0, const double eps0)
    {
        double rxy = sqrt(x.x*x.x+x.y*x.y);
        double tmp = rxy-R0;
        if (tmp*tmp + x.z*x.z - eps0*eps0 < 0)
        {
            double f = strength/(4*M_PI*R0) * (log(8.*R0/eps0) -.25 + eps0/R0);
            return Point3d(-x.y/rxy*f, x.x/rxy*f, 0 );
            //return Point3d(1,1,1);
        }
        else
        {
            return Point3d(0);
        }
    }


    Point3d omega_DoubleVortexRing_z(const Point3d x, const double strength)
    {
        double dist = 2.0;
        Point3d c_a = Point3d(0.,0.,-dist/2.);
        Point3d c_b = Point3d(0.,0.,+dist/2.);
        return omega_VortexRing_z(x-c_a, strength) - omega_VortexRing_z(x-c_b, strength);
    }

    Point3d omega_LambOseen_z(const Point3d x, const double strength, const double core_radius)
    {
        double r2 = x.x*x.x+x.y*x.y;//dot(x,x);
        double rct_squared = 0.15*0.15;
        return Point3d(0,0,strength/(M_PI*rct_squared) * exp(-r2/rct_squared));
    }
    Point3d omega_CounterRotatingVortexPair_z(const Point3d x, const double dist, const double strength, const double core_radius)
    {
        Point3d c_a = Point3d(-dist/2.,0.,0.);
        Point3d c_b = Point3d(+dist/2.,0.,0.);
        return - omega_LambOseen_z(x-c_a, strength, core_radius) + omega_LambOseen_z(x-c_b, strength, core_radius);
    }
    Point3d omega_CoRotatingVortexPair_z(const Point3d x, const double dist, const double strength, const double core_radius)
    {
        Point3d c_a = Point3d(-dist/2.,0.,0.);
        Point3d c_b = Point3d(+dist/2.,0.,0.);
        return omega_LambOseen_z(x-c_a, strength, core_radius) + omega_LambOseen_z(x-c_b, strength, core_radius);
    }

    Point3d omega_LambOseen_x(const Point3d x, const double strength, const double core_radius)
    {
        double r2 = x.y*x.y+x.z*x.z;//dot(x,x);
        double rct_squared = 0.15*0.15;
        return Point3d(strength/(M_PI*rct_squared) * exp(-r2/rct_squared),0,0);
    }

    Point3d omega_CounterRotatingVortexPair_x(const Point3d x, const double dist, const double strength, const double core_radius)
    {
        Point3d c_a = Point3d(0.,-dist/2.,0.);
        Point3d c_b = Point3d(0.,+dist/2.,0.);
        return - omega_LambOseen_x(x-c_a, strength, core_radius) + omega_LambOseen_x(x-c_b, strength, core_radius);
    }
    Point3d omega_CoRotatingVortexPair_x(const Point3d x, const double dist, const double strength, const double core_radius)
    {
        Point3d c_a = Point3d(0.,-dist/2.,0.);
        Point3d c_b = Point3d(0.,+dist/2.,0.);
        return omega_LambOseen_x(x-c_a, strength, core_radius) + omega_LambOseen_x(x-c_b, strength, core_radius);
    }


    Point3d omega_LambOseen_y(const Point3d x, const double strength, const double core_radius)
    {
        double r2 = x.x*x.x+x.z*x.z;//dot(x,x);
        double rct_squared = 0.15*0.15;
        return Point3d(0,strength/(M_PI*rct_squared) * exp(-r2/rct_squared),0);
    }
    Point3d omega_CounterRotatingVortexPair_y(const Point3d x, const double dist, const double strength, const double core_radius)
    {
        Point3d c_a = Point3d(0.,0.,-dist/2.);
        Point3d c_b = Point3d(0.,0.,+dist/2.);
        return - omega_LambOseen_y(x-c_a, strength, core_radius) + omega_LambOseen_y(x-c_b, strength, core_radius);
    }
    Point3d omega_CoRotatingVortexPair_y(const Point3d x, const double dist, const double strength, const double core_radius)
    {
        Point3d c_a = Point3d(0.,0.,-dist/2.);
        Point3d c_b = Point3d(0.,0.,+dist/2.);
        return omega_LambOseen_y(x-c_a, strength, core_radius) + omega_LambOseen_y(x-c_b, strength, core_radius);
    }

    void init(
            const std::shared_ptr<VPM::Parameters> params,
            const int example_num,
            const double dist,
            const double strength,
            const double core_radius,
            std::vector<Point3d> & positions,
            std::vector<Point3d> & omega
            )
    {

        positions.clear();
        omega.clear();

        int n = -1;
        for (unsigned int npz=0; npz<params->m_num_pz; npz++)
        {
            for (unsigned int npy=0; npy<params->m_num_py; npy++)
            {
                for (unsigned int npx=0; npx<params->m_num_px; npx++)
                {
                    Point3d pos = params->m_domain_ll + Point3d(
                            double(2*npx+1)/double(2*params->m_num_px)*(params->m_domain_ur.x-params->m_domain_ll.x),
                            double(2*npy+1)/double(2*params->m_num_py)*(params->m_domain_ur.y-params->m_domain_ll.y),
                            double(2*npz+1)/double(2*params->m_num_pz)*(params->m_domain_ur.z-params->m_domain_ll.z)
                            );
                    Point3d w;
                    switch (example_num)
                    {
                        case -1:
                            w = Point3d(0.);
                            break;
                        case 0:
                            w = omega_CounterRotatingVortexPair_x(pos, dist, strength, core_radius);
                            break;
                        case 1:
                            w = omega_CounterRotatingVortexPair_y(pos, dist, strength, core_radius);
                            break;
                        case 2:
                            w = omega_CounterRotatingVortexPair_z(pos, dist, strength, core_radius);
                            break;
                        case 3:
                            w = omega_CoRotatingVortexPair_x(pos, dist, strength, core_radius);
                            break;
                        case 4:
                            w = omega_CoRotatingVortexPair_y(pos, dist, strength, core_radius);
                            break;
                        case 5:
                            w = omega_CoRotatingVortexPair_z(pos, dist, strength, core_radius);
                            break;
                        case 6:
                            {
                                double dist = 2.0;
                                Point3d c_a = Point3d(0.,0.,-dist/2.);
                                w = omega_VortexRing_z(pos-c_a,strength);
                            }
                            break;
                        case 7:
                            w = omega_DoubleVortexRing_z(pos,strength);
                            break;
                    }
                    positions.push_back(pos);
                    omega.push_back(w);
                }
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
            std::cerr<<"---> Writing particles: Warning: Velocity field does not correspond to vorticity"<<std::endl;
        }
        if (!pf.cartesianGrid)
        {
            std::cerr<<"---> Writing particles: Warning: Not on Cartesian grid"<<std::endl;
        }

        std::string t = std::to_string(timestep);

        pf.params->writeToFile(filename+"_t"+t+"_conf.dat");

        //std::ofstream myfile_pos;
        std::ofstream myfile_omega;
        std::ofstream myfile_vel;
        std::ofstream myfile_Jvel_x_omaga;

        //myfile_pos.open (filename+"_t"+t+"_pos.dat", std::ofstream::binary);
        myfile_omega.open (filename+"_t"+t+"_omega.dat", std::ofstream::binary);
        myfile_vel.open (filename+"_t"+t+"_U.dat", std::ofstream::binary);
        myfile_Jvel_x_omaga.open (filename+"_t"+t+"_Jxw.dat", std::ofstream::binary);

        for (unsigned int i=0; i<pf.positions.size(); i++)
        {
            //myfile_pos.write(reinterpret_cast<const char*>(&pf.positions[i].x), sizeof(double));
            //myfile_pos.write(reinterpret_cast<const char*>(&pf.positions[i].y), sizeof(double));
            //myfile_pos.write(reinterpret_cast<const char*>(&pf.positions[i].z), sizeof(double));
            myfile_omega.write(reinterpret_cast<const char*>(&pf.omega[i].x), sizeof(double));
            myfile_omega.write(reinterpret_cast<const char*>(&pf.omega[i].y), sizeof(double));
            myfile_omega.write(reinterpret_cast<const char*>(&pf.omega[i].z), sizeof(double));
            myfile_vel.write(reinterpret_cast<const char*>(&pf.velocity[i].x), sizeof(double));
            myfile_vel.write(reinterpret_cast<const char*>(&pf.velocity[i].y), sizeof(double));
            myfile_vel.write(reinterpret_cast<const char*>(&pf.velocity[i].z), sizeof(double));
            myfile_Jvel_x_omaga.write(reinterpret_cast<const char*>(&pf.Jvel_x_omega[i].x), sizeof(double));
            myfile_Jvel_x_omaga.write(reinterpret_cast<const char*>(&pf.Jvel_x_omega[i].y), sizeof(double));
            myfile_Jvel_x_omaga.write(reinterpret_cast<const char*>(&pf.Jvel_x_omega[i].z), sizeof(double));
        }

        //myfile_pos.close();
        myfile_omega.close();
        myfile_vel.close();
        myfile_Jvel_x_omaga.close();

    }

    bool readParticlesFromFile(
            const std::string & filename,
            const std::shared_ptr<VPM::Parameters> params,
            std::vector<Point3d> & positions,
            std::vector<Point3d> & omega
            )
    {
        //std::ifstream myfile_pos(filename+"_pos.dat", std::ifstream::binary);
        //if(!myfile_pos.is_open())
        //{
        //    std::cerr<<"\nCannot open file: "<< filename<<"_pos.dat";
        //    return false;
        //}
        std::ifstream myfile_omega(filename+"_omega.dat", std::ifstream::binary);
        if(!myfile_omega.is_open())
        {
            std::cerr<<"\nCannot open file: "<< filename<<"_omega.dat";
            return false;
        }

        omega.clear();

        for (unsigned int i=0; i<params->m_N; i++)
        {
            double tmp_x, tmp_y, tmp_z;

            //myfile_pos.read((char *) &tmp_x, sizeof(double));
            //myfile_pos.read((char *) &tmp_y, sizeof(double));
            //myfile_pos.read((char *) &tmp_z, sizeof(double));
            //positions.push_back(Point3d(tmp_x, tmp_y, tmp_y));

            myfile_omega.read((char *) &tmp_x, sizeof(double));
            myfile_omega.read((char *) &tmp_y, sizeof(double));
            myfile_omega.read((char *) &tmp_z, sizeof(double));
            omega.push_back(Point3d(tmp_x, tmp_y, tmp_y));
        }

        myfile_omega.close();

        positions.clear();
        for (unsigned int npz=0; npz<params->m_num_pz; npz++)
        {
            for (unsigned int npy=0; npy<params->m_num_py; npy++)
            {
                for (unsigned int npx=0; npx<params->m_num_px; npx++)
                {
                    Point3d pos = params->m_domain_ll + Point3d(
                            double(2*npx+1)/double(2*params->m_num_px)*(params->m_domain_ur.x-params->m_domain_ll.x),
                            double(2*npy+1)/double(2*params->m_num_py)*(params->m_domain_ur.y-params->m_domain_ll.y),
                            double(2*npz+1)/double(2*params->m_num_pz)*(params->m_domain_ur.z-params->m_domain_ll.z)
                            );
                    positions.push_back(pos);
                }
            }
        }

        return true;
    }

}
