#pragma once
#include "Point2d.hpp"

#include <iostream>
#include <fstream>
#include <memory>

#define PSEorder 2

namespace VPM
{

    enum interpolationMethod { Gamma3=1, M4d=2, M5d=3, M6d=4, /*order6,*/ Gamma3_fast=5, M4d_fast=6, M5d_fast=7, M6d_fast=8/*, order6_fast*/};

    struct RemeshParams
    {
        bool m_isOn;
        int m_steps;
        std::string m_methodstring;
        int m_method;

        RemeshParams(){};

        RemeshParams(
                const bool isOn = true,
                const int steps = 1,
                const std::string methodstring = "M4d_fast"
              )
            :
                m_isOn(isOn),
                m_steps(steps),
                m_methodstring(methodstring)
        {
            if (m_methodstring=="Gamma3")
            {
                m_method = Gamma3;
            }
            else if (m_methodstring=="M4d")
            {
                m_method = M4d;
            }
            else if (m_methodstring=="M5d")
            {
                m_method = M5d;
            }
            else if (m_methodstring=="M6d")
            {
                m_method = M6d;
            }
            //else if (m_methodstring=="order6")
            //{
            //    m_method = order6;
            //}
            else if (m_methodstring=="Gamma3_fast")
            {
                m_method = Gamma3_fast;
            }
            else if (m_methodstring=="M4d_fast")
            {
                m_method = M4d_fast;
            }
            else if (m_methodstring=="M5d_fast")
            {
                m_method = M5d_fast;
            }
            else if (m_methodstring=="M6d_fast")
            {
                m_method = M6d_fast;
            }
            //else if (m_methodstring=="order6_fast")
            //{
            //    m_method = order6_fast;
            //}
            else
            {
                std::cerr<<"\n ---> ERROR!!! Interpolation method \""<<m_methodstring<<"\" is not yet implemented!";
            }
        }

        RemeshParams(
                const bool isOn = true,
                const int steps = 1,
                const int method = M4d_fast
              )
            :
                m_isOn(isOn),
                m_steps(steps),
                m_method(method)
        {
            switch (m_method)
            {
                case(Gamma3):
                    m_methodstring = "Gamma3";
                    break;
                case(M4d):
                    m_methodstring = "M4d";
                    break;
                case(M5d):
                    m_methodstring = "M5d";
                    break;
                case(M6d):
                    m_methodstring = "M6d";
                    break;
                //case(order6):
                //    m_methodstring = "order6";
                //    break;
                case(Gamma3_fast):
                    m_methodstring = "Gamma3_fast";
                    break;
                case(M4d_fast):
                    m_methodstring = "M4d_fast";
                    break;
                case(M5d_fast):
                    m_methodstring = "M5d_fast";
                    break;
                case(M6d_fast):
                    m_methodstring = "M6d_fast";
                    break;
                //case(order6_fast):
                //    m_methodstring = "order6_fast";
                //    break;
                default:
                    std::cerr<<"\n ---> ERROR!!! Interpolation method \""<<m_methodstring<<"\" is not yet implemented!";
                    break;
            }
        }

    };


    enum boundary_conditions {
        BC_XL_NONE, BC_XL_WALL, BC_XL_OUTFLOW,   BC_XR_NONE, BC_XR_WALL, BC_XR_OUTFLOW,
        BC_YL_NONE, BC_YL_WALL, BC_YL_OUTFLOW,   BC_YR_NONE, BC_YR_WALL, BC_YR_OUTFLOW,
        BC_ZL_NONE, BC_ZL_WALL, BC_ZL_OUTFLOW,   BC_ZR_NONE, BC_ZR_WALL, BC_ZR_OUTFLOW
    };

    struct BCParams
    {
        int m_XL, m_XR;
        int m_YL, m_YR;
        int m_ZL, m_ZR;
        double m_threshold_xl, m_threshold_xr;
        double m_threshold_yl, m_threshold_yr;
        double m_threshold_zl, m_threshold_zr;

        BCParams(){};

        BCParams(
                const int xl = BC_XL_NONE,
                const double to_xl = double(0),
                const int xr = BC_XR_NONE,
                const double to_xr = double(0),
                const int yl = BC_YL_NONE,
                const double to_yl = double(0),
                const int yr = BC_YR_NONE,
                const double to_yr = double(0),
                const int zl = BC_ZL_NONE,
                const double to_zl = double(0),
                const int zr = BC_ZR_NONE,
                const double to_zr = double(0)
              )
            :
                m_XL(xl),
                m_XR(xr),
                m_YL(yl),
                m_YR(yr),
                m_ZL(zl),
                m_ZR(zr),
                m_threshold_xl(to_xl),
                m_threshold_xr(to_xr),
                m_threshold_yl(to_yl),
                m_threshold_yr(to_yr),
                m_threshold_zl(to_zl),
                m_threshold_zr(to_zr)
        {
        }

        BCParams(
                const std::string xlstring = "none",
                const double to_xl = double(0),
                const std::string xrstring = "none",
                const double to_xr = double(0),
                const std::string ylstring = "none",
                const double to_yl = double(0),
                const std::string yrstring = "none",
                const double to_yr = double(0),
                const std::string zlstring = "none",
                const double to_zl = double(0),
                const std::string zrstring = "none",
                const double to_zr = double(0)
              )
            :
                m_threshold_xl(to_xl),
                m_threshold_xr(to_xr),
                m_threshold_yl(to_yl),
                m_threshold_yr(to_yr),
                m_threshold_zl(to_zl),
                m_threshold_zr(to_zr)
        {
            // XL
            if (xlstring == "none")
            {
                m_XL = BC_XL_NONE;
            }
            else if (xlstring== "wall")
            {
                m_XL = BC_XL_WALL;
            }
            else if (xlstring == "outflow")
            {
                m_XL = BC_XL_OUTFLOW;
            }
            else
            {
                std::cerr<<"\n ---> ERROR!!! BC \""<<xlstring<<"\" is not yet implemented!";
            }
            // XR
            if (xrstring == "none")
            {
                m_XR = BC_XR_NONE;
            }
            else if (xrstring== "wall")
            {
                m_XR = BC_XR_WALL;
            }
            else if (xrstring == "outflow")
            {
                m_XR = BC_XR_OUTFLOW;
            }
            else
            {
                std::cerr<<"\n ---> ERROR!!! BC \""<<xrstring<<"\" is not yet implemented!";
            }
            // YL
            if (ylstring == "none")
            {
                m_YL = BC_YL_NONE;
            }
            else if (ylstring== "wall")
            {
                m_YL = BC_YL_WALL;
            }
            else if (ylstring == "outflow")
            {
                m_YL = BC_YL_OUTFLOW;
            }
            else
            {
                std::cerr<<"\n ---> ERROR!!! BC \""<<ylstring<<"\" is not yet implemented!";
            }
            // YR
            if (yrstring == "none")
            {
                m_YR = BC_YR_NONE;
            }
            else if (yrstring== "wall")
            {
                m_YR = BC_YR_WALL;
            }
            else if (yrstring == "outflow")
            {
                m_YR = BC_YR_OUTFLOW;
            }
            else
            {
                std::cerr<<"\n ---> ERROR!!! BC \""<<yrstring<<"\" is not yet implemented!";
            }
            // ZL
            if (zlstring == "none")
            {
                m_ZL = BC_ZL_NONE;
            }
            else if (zlstring== "wall")
            {
                m_ZL = BC_ZL_WALL;
            }
            else if (zlstring == "outflow")
            {
                m_ZL = BC_ZL_OUTFLOW;
            }
            else
            {
                std::cerr<<"\n ---> ERROR!!! BC \""<<zlstring<<"\" is not yet implemented!";
            }
            // ZR
            if (zrstring == "none")
            {
                m_ZR = BC_ZR_NONE;
            }
            else if (zrstring== "wall")
            {
                m_ZR = BC_ZR_WALL;
            }
            else if (zrstring == "outflow")
            {
                m_ZR = BC_ZR_OUTFLOW;
            }
            else
            {
                std::cerr<<"\n ---> ERROR!!! BC \""<<zrstring<<"\" is not yet implemented!";
            }
        }

    };
    struct Parameters
    {

        double m_nu;
        double m_nu_lambda_bydt;
        double m_nu_lambda_sq_bydt;
        unsigned int m_num_px;
        unsigned int m_regular_num_px;
        unsigned int m_num_py;
        unsigned int m_regular_num_py;
        VPM::Point2d m_domain_ll;
        VPM::Point2d m_domain_ur;

        VPM::Point2d m_Uinfty;

        VPM::Point2d m_initial_domain_ll;
        VPM::Point2d m_initial_domain_ur;
        VPM::Point2d m_shift_ll;
        VPM::Point2d m_scale;

        std::shared_ptr<RemeshParams> m_remesh;
        std::shared_ptr<BCParams> m_bc;

        unsigned int m_order_ODEsolver;

        double m_population_threshold;
        double m_population_threshold2;

        double m_dx;
        double m_dy;
        double m_vol;
        double m_eps;
        double m_sigma;
        double m_sigma2;
        unsigned int m_N;// Number of charges;

        bool m_isInitialized;

        Parameters()
        {
        };

        Parameters(
                const VPM::Point2d domain_ll,
                const VPM::Point2d domain_ur,
                const unsigned int num_px,
                const unsigned int num_py,
                const double nu,
                const double population_threshold,
                const std::shared_ptr<RemeshParams> remesh,
                const std::shared_ptr<BCParams> bc,
                const unsigned int order_ODEsolver,
                const VPM::Point2d Uinfty
                )
            :
                m_domain_ll(domain_ll),
                m_domain_ur(domain_ur),
                m_initial_domain_ll(domain_ll),
                m_initial_domain_ur(domain_ur),
                m_num_px(num_px),
                m_regular_num_px(num_px),
                m_num_py(num_py),
                m_regular_num_py(num_py),
                m_nu(nu),
                m_population_threshold(population_threshold),
                m_order_ODEsolver(order_ODEsolver),
                m_Uinfty(Uinfty)
        {
            m_remesh = remesh;
            m_bc = bc;
            m_dx = (m_domain_ur.x-m_domain_ll.x)/double(m_num_px);
            m_dy = (m_domain_ur.y-m_domain_ll.y)/double(m_num_py);
            m_eps = 2*std::min(m_dx, m_dy);
            m_sigma = PSEorder*std::min(m_dx, m_dy);
            m_sigma2 = m_sigma*m_sigma;
            m_population_threshold2 = m_population_threshold*m_population_threshold;
            m_vol = m_dx*m_dy;
            m_N = m_num_px*m_num_py;

            if (abs(m_dx-m_dy)>1e-6)
            {
                std::cerr<<"\nWarning!!! Implementation assumes (at the moment dx=dy)!!!, but: dx="<<m_dx<<"!="<<m_dy<<"=dy and px="<<m_num_px<<", py="<<m_num_py;
                std::cerr<<"\nWarning!!! Implementation assumes (at the moment dx=dy)!!!, but: dx="<<m_dx<<"!="<<m_dy<<"=dy and px="<<m_num_px<<", py="<<m_num_py;
                std::cerr<<"\nWarning!!! Implementation assumes (at the moment dx=dy)!!!, but: dx="<<m_dx<<"!="<<m_dy<<"=dy and px="<<m_num_px<<", py="<<m_num_py;
            }
            m_scale = VPM::Point2d(1,1);
            m_shift_ll = VPM::Point2d(0,0);
            m_isInitialized = true;
        };

        void print() const
        {
            std::cerr<<"nu=" << std::to_string(m_nu) <<std::endl;
            std::cerr<<"nulambda=" << std::to_string(m_nu_lambda_bydt) <<std::endl;
            std::cerr<<"nulambda_sq=" << std::to_string(m_nu_lambda_sq_bydt) <<std::endl;
            std::cerr<<"nx=" << std::to_string(m_num_px) <<std::endl;
            std::cerr<<"regnx=" << std::to_string(m_regular_num_px) <<std::endl;
            std::cerr<<"ny=" << std::to_string(m_num_py) <<std::endl;
            std::cerr<<"regny=" << std::to_string(m_regular_num_py) <<std::endl;
            std::cerr<<"domainll=" << std::to_string(m_domain_ll.x) <<","<< std::to_string(m_domain_ll.y) <<std::endl;
            std::cerr<<"domainur=" << std::to_string(m_domain_ur.x) <<","<< std::to_string(m_domain_ur.y) <<std::endl;
            std::cerr<<"initialdomainll=" << std::to_string(m_initial_domain_ll.x) <<","<< std::to_string(m_initial_domain_ll.y) <<std::endl;
            std::cerr<<"initialdomainur=" << std::to_string(m_initial_domain_ur.x) <<","<< std::to_string(m_initial_domain_ur.y) <<std::endl;
            std::cerr<<"Uinfty=" << std::to_string(m_Uinfty.x) <<","<< std::to_string(m_Uinfty.y) <<std::endl;
            std::cerr<<"shift_ll=" << std::to_string(m_shift_ll.x) <<","<< std::to_string(m_shift_ll.y) <<std::endl;
            std::cerr<<"scale=" << std::to_string(m_scale.x) <<","<< std::to_string(m_scale.y) <<std::endl;
//
//            std::cerr<<"" << std::to_string(std::shared_ptr<RemeshParams> m_remesh) <<std::endl;
//            std::cerr<<"" << std::to_string(std::shared_ptr<BCParams> m_bc) <<std::endl;
//
            std::cerr<<"odesolver=" << std::to_string(m_order_ODEsolver) <<std::endl;

            std::cerr<<"pop_threshold=" << std::to_string(m_population_threshold) <<std::endl;
            std::cerr<<"pop_threshold2=" << std::to_string(m_population_threshold2) <<std::endl;

            std::cerr<<"dx=" << std::to_string(m_dx) <<std::endl;
            std::cerr<<"dy=" << std::to_string(m_dy) <<std::endl;
            std::cerr<<"vol=" << std::to_string(m_vol) <<std::endl;
            std::cerr<<"eps=" << std::to_string(m_eps) <<std::endl;
            std::cerr<<"sigma=" << std::to_string(m_sigma) <<std::endl;
            std::cerr<<"sigma2=" << std::to_string(m_sigma2) <<std::endl;
            std::cerr<<"N=" << std::to_string(m_N) <<std::endl;// Number of charges;

            std::cerr<<"isinitialized=" << std::to_string(m_isInitialized) <<std::endl;
        }

        void updateDomain(const VPM::Point2d & new_domain_ll, const VPM::Point2d & new_domain_ur)
        {
            m_scale = new_domain_ur-new_domain_ll;
            m_scale.x /= m_initial_domain_ur.x-m_initial_domain_ll.x;
            m_scale.y /= m_initial_domain_ur.y-m_initial_domain_ll.y;

            //std::cerr<<"\n initialdomain="<<m_initial_domain_ll.x<<", "<<m_initial_domain_ur.x;
            //std::cerr<<"\n newdomain X ="<<new_domain_ll.x<<", "<<new_domain_ur.x;

            m_shift_ll = new_domain_ll-m_initial_domain_ll;
            m_shift_ll.x /= m_initial_domain_ur.x-m_initial_domain_ll.x;
            m_shift_ll.y /= m_initial_domain_ur.y-m_initial_domain_ll.y;

            //std::cerr<<"\n scale="<<m_scale.x<<", "<<m_scale.y<<", shift="<<m_shift_ll.x<<", "<<m_shift_ll.y;

            m_domain_ll = new_domain_ll;
            m_domain_ur = new_domain_ur;
        }


        bool readFromFile(std::string filename)
        {
            std::ifstream myfile_conf(filename, std::ifstream::binary);
            if(!myfile_conf.is_open())
            {
                std::cerr<<"\nCannot open file: "<< filename;
                m_isInitialized = false;
                return false;
            }

            int remesh_isOn;
            int remesh_steps;
            unsigned int remesh_method;

            int bc_XL;
            int bc_XR;
            int bc_YL;
            int bc_YR;
            double bc_to_XL;
            double bc_to_XR;
            double bc_to_YL;
            double bc_to_YR;

            myfile_conf.read((char *) &m_domain_ll.x, sizeof(double));
            myfile_conf.read((char *) &m_domain_ll.y, sizeof(double));
            myfile_conf.read((char *) &m_domain_ur.x, sizeof(double));
            myfile_conf.read((char *) &m_domain_ur.y, sizeof(double));
            //myfile_conf.read((char *) &eps, sizeof(double));
            myfile_conf.read((char *) &m_nu, sizeof(double));
            //myfile_conf.read((char *) &m_final_time, sizeof(double));
            //myfile_conf.read((char *) &sigma, sizeof(double));
            myfile_conf.read((char *) &m_population_threshold, sizeof(double));
            myfile_conf.read((char *) &bc_to_XL, sizeof(double));
            myfile_conf.read((char *) &bc_to_XR, sizeof(double));
            myfile_conf.read((char *) &bc_to_YL, sizeof(double));
            myfile_conf.read((char *) &bc_to_YR, sizeof(double));
            myfile_conf.read((char *) &m_Uinfty.x, sizeof(double));
            myfile_conf.read((char *) &m_Uinfty.y, sizeof(double));

            // next 3 are remesh parameters
            myfile_conf.read((char *) &remesh_isOn, sizeof(int));
            myfile_conf.read((char *) &remesh_steps, sizeof(int));
            myfile_conf.read((char *) &remesh_method, sizeof(int));
            //
            myfile_conf.read((char *) &m_N, sizeof(unsigned int));
            myfile_conf.read((char *) &m_num_px, sizeof(unsigned int));
            myfile_conf.read((char *) &m_num_py, sizeof(unsigned int));
            myfile_conf.read((char *) &m_order_ODEsolver, sizeof(unsigned int));

            myfile_conf.read((char *) &bc_XL, sizeof(int));
            myfile_conf.read((char *) &bc_XR, sizeof(int));
            myfile_conf.read((char *) &bc_YL, sizeof(int));
            myfile_conf.read((char *) &bc_YR, sizeof(int));

            myfile_conf.close();


            m_regular_num_px = m_num_px;
            m_regular_num_py = m_num_py;
            m_initial_domain_ll = m_domain_ll;
            m_initial_domain_ur = m_domain_ur;

            m_remesh = std::make_shared<RemeshParams>
                (bool(remesh_isOn), remesh_steps, remesh_method);

            m_bc = std::make_shared<BCParams>
                (bc_XL, bc_to_XL, bc_XR, bc_to_XR, bc_YL, bc_to_YL, bc_YR, bc_to_YR
                 );

            m_dx = (m_domain_ur.x-m_domain_ll.x)/double(m_num_px);
            m_dy = (m_domain_ur.y-m_domain_ll.y)/double(m_num_py);
            m_eps = 2*std::min(m_dx, m_dy);
            m_sigma = PSEorder*std::min(m_dx, m_dy);
            m_sigma2 = m_sigma*m_sigma;
            m_population_threshold2 = m_population_threshold*m_population_threshold;
            m_vol = m_dx*m_dy;
            m_N = m_num_px*m_num_py;

            if (abs(m_dx-m_dy)>1e-6)
            {
                std::cerr<<"\nWarning!!! Implementation assumes (at the moment dx=dy)!!!, but: dx="<<m_dx<<"!="<<m_dy<<"=dy and px="<<m_num_px<<", py="<<m_num_py;
                std::cerr<<"\nWarning!!! Implementation assumes (at the moment dx=dy)!!!, but: dx="<<m_dx<<"!="<<m_dy<<"=dy and px="<<m_num_px<<", py="<<m_num_py;
                std::cerr<<"\nWarning!!! Implementation assumes (at the moment dx=dy)!!!, but: dx="<<m_dx<<"!="<<m_dy<<"=dy and px="<<m_num_px<<", py="<<m_num_py;
            }
            m_scale = VPM::Point2d(1,1);
            m_shift_ll = VPM::Point2d(0,0);
            m_isInitialized = true;

            return true;
        };

        void writeToFile(std::string filename) const
        {

            std::ofstream myfile_conf(filename, std::ofstream::binary);

            myfile_conf.write((char *) &m_domain_ll.x, sizeof(double));
            myfile_conf.write((char *) &m_domain_ll.y, sizeof(double));
            myfile_conf.write((char *) &m_domain_ur.x, sizeof(double));
            myfile_conf.write((char *) &m_domain_ur.y, sizeof(double));
            //myfile_conf.write((char *) &m_eps, sizeof(double));
            myfile_conf.write((char *) &m_nu, sizeof(double));
            //myfile_conf.write((char *) &m_final_time, sizeof(double));
            //myfile_conf.write((char *) &m_sigma, sizeof(double));
            myfile_conf.write((char *) &m_population_threshold, sizeof(double));
            myfile_conf.write((char *) &m_bc->m_threshold_xl, sizeof(double));
            myfile_conf.write((char *) &m_bc->m_threshold_xr, sizeof(double));
            myfile_conf.write((char *) &m_bc->m_threshold_yl, sizeof(double));
            myfile_conf.write((char *) &m_bc->m_threshold_yr, sizeof(double));
            myfile_conf.write((char *) &m_Uinfty.x, sizeof(double));
            myfile_conf.write((char *) &m_Uinfty.y, sizeof(double));

            // next 3 are remesh parameters
            int remeshison = int(m_remesh->m_isOn);
            myfile_conf.write((char *) &remeshison, sizeof(int));
            myfile_conf.write((char *) &m_remesh->m_steps, sizeof(int));
            myfile_conf.write((char *) &m_remesh->m_method, sizeof(int));
            //
            myfile_conf.write((char *) &m_N, sizeof(unsigned int));
            myfile_conf.write((char *) &m_num_px, sizeof(unsigned int));
            myfile_conf.write((char *) &m_num_py, sizeof(unsigned int));
            myfile_conf.write((char *) &m_order_ODEsolver, sizeof(unsigned int));
            myfile_conf.write((char *) &m_bc->m_XL, sizeof(int));
            myfile_conf.write((char *) &m_bc->m_XR, sizeof(int));
            myfile_conf.write((char *) &m_bc->m_YL, sizeof(int));
            myfile_conf.write((char *) &m_bc->m_YR, sizeof(int));

            myfile_conf.close();
        };


    };

}
