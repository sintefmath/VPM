#include "Redistribute.hpp"
#include "FlattenArrayIndex.hpp"

#include <math.h>

namespace VPM
{

    namespace {


        inline double Gamma3_dir(double p)
        {
            p = std::abs(p);
            if (p>=2)
            {
                return 0;
            }
            else if (p>=1)
            {
                return 1./6.*(1-p)*(2-p)*(3-p);
            }
            else // p<1
            {
                return 0.5*(1-p*p)*(2-p);
            }

        }
        inline double interpolation_Kernel_Gamma3(const VPM::Point2d & p, const double delta)
        {
            VPM::Point2d pd = p/delta;
            return Gamma3_dir(pd.x)*Gamma3_dir(pd.y);
        }

        inline double M4d_dir(double p)
        {
            p = std::abs(p);
            if (p>=2)
            {
                return 0;
            }
            else if (p>=1)
            {
                double twomp = 2.-p;
                double onemp = 1.-p;
                return 0.5*twomp*twomp*onemp;
            }
            else // p<1
            {
                double p2 = p*p;
                return 1.-5./2.*p2+3./2.*p2*p;
            }

        }
        inline double interpolation_Kernel_M4d(const VPM::Point2d & p, const double delta)
        {
            VPM::Point2d pd = p;
            return M4d_dir(pd.x/delta)*M4d_dir(pd.y/delta);
        }

        inline double M5d_dir(double p)
        {
            p = std::abs(p);
            if (p>=2)
            {
                return 0;
            }
            else if (p>=1)
            {
                double p2 = p*p;
                return 1. - 11./6.*p + p2 - p2*p/6.;
            }
            else // p<1
            {
                double p2 = p*p;
                return 1. - 0.5*p - p2 + 0.5*p2*p;
            }

        }
        inline double interpolation_Kernel_M5d(const VPM::Point2d & p, const double delta)
        {
            VPM::Point2d pd = p/delta;
            return M5d_dir(pd.x)*M5d_dir(pd.y);
        }

        inline double M6d_dir(double p)
        {
            p = std::abs(p);
            if (p>=3)
            {
                return 0;
            }
            else if (p>=2)
            {
                double p2 = p*p;
                double p4 = p2*p2;
                return 1. - 137./60.*p + 15./8.*p2 - 17./24.*p2*p + 1./8.*p4 - 1./120.*p4*p;
            }
            else if (p>=1)
            {
                double p2 = p*p;
                double p4 = p2*p2;
                return 1. - 13./12.*p - 5./8.*p2 + 25./24.*p2*p - 3./8.*p4 + 1./24*p4*p;
            }
            else // p<1
            {
                double p2 = p*p;
                double p4 = p2*p2;
                return 1. - 1./3.*p - 5./4.*p2 + 5./12.*p2*p + 1./4.*p4 - 1./12.*p4*p;
            }

        }
        inline double interpolation_Kernel_M6d(const VPM::Point2d & p, const double delta)
        {
            VPM::Point2d pd = p/delta;
            return M6d_dir(pd.x)*M6d_dir(pd.y);
        }

        inline double interpolation_Kernel_order6(const VPM::Point2d & p, const double delta)
        {
            VPM::Point2d pd = p/delta;

            double px2 = pd.x*pd.x;
            double py2 = pd.y*pd.y;
            double r2 = px2 + py2;
            if (r2>=16)
            {
                return 0;
            }

            double tmp_x = 15./4. - 5.*px2 + px2*px2;
            double tmp_y = 15./4. - 5.*py2 + py2*py2;
            double res = 1./M_PI*tmp_x/2.*tmp_y/2.*exp(-r2);

            return res;
        }

    }

    Redistribute::Redistribute()
    {
    }
    Redistribute::~Redistribute()
    {
    }

    void Redistribute::get_indices(
            const Parameters& params,
            const int stenc,
            const VPM::Point2d & pos,
            VPM::IPoint2d & ind_min,
            VPM::IPoint2d & ind_max
            )
    {

        VPM::Point2d scaled_pos = (pos - params.m_domain_ll)/(params.m_domain_ur - params.m_domain_ll);
        int i_ind = std::floor(scaled_pos.x*double(params.m_num_px-1));
        ind_min.x = std::max(int(0), i_ind-stenc);
        ind_max.x = std::min(int(params.m_num_px-1), i_ind+stenc+1);
        int j_ind = std::floor(scaled_pos.y*double(params.m_num_py-1));
        ind_min.y = std::max(int(0), j_ind-stenc);
        ind_max.y = std::min(int(params.m_num_py-1), j_ind+stenc+1);
    }

    void Redistribute::redistribute(
            ParticleField & pf,
            const REDISTRIBUTE2D g2p
            )
    {
        redistribute(
                pf.positions,
                pf.regular_positions,
                pf.params,
                pf.omega,
                g2p
                );
        pf.cartesianGrid = true;
    }

    template <class T>
    void Redistribute::redistribute(
            const std::vector<Point2d> & old_pos,
            const std::vector<Point2d> & new_pos,
            Parameters & params,
            std::vector<T> & field,
            const REDISTRIBUTE2D g2p
            )
    {
        // 1) copy field to field_old
        std::vector<T> field_old = field;
        // 2) initialize field to 0
        unsigned long new_N = new_pos.size();
        field.clear();
        field.resize(new_N, T(0.));

        if (g2p == RED2D_G2P_BUILD)
        {
            m_grid_to_particles.clear();
            m_grid_to_particles.resize(new_N);
        }

        // declaring this const, helps the compiler to optimize (decision outside the for loop(s))
        const auto method = params.m_remesh->m_method;

        std::string method_str;
        switch (method)
        {
            case(VPM::Gamma3):
                method_str = "Gamma3";
                break;
            case(VPM::M4d):
                method_str = "M4d";
                break;
            case(VPM::M5d):
                method_str = "M5d";
                break;
            case(VPM::M6d):
                method_str = "M6d";
                break;
            case(VPM::Gamma3_fast) :
                method_str = "Gamma3_fast";
                break;
            case(VPM::M4d_fast) :
                method_str = "M4d_fast";
                break;
            case(VPM::M5d_fast) :
                method_str = "M5d_fast";
                break;
            case(VPM::M6d_fast):
                method_str = "M6d_fast";
                break;
            default:
                std::cerr<<"    ---> remesh method not found!!!\n";
                exit(0);
        }
        #ifdef VPM_VERBOSE
        std::cerr<<"    ---> redistribution method: "<<method_str <<std::endl;
        #endif
        switch (method)
        {
            case(VPM::Gamma3): case(VPM::M4d): case(VPM::M5d): case(VPM::M6d):
                {
                    for(unsigned int j=0; j<new_N; j++ )
                    {
                        for(unsigned int i=0; i<params.m_N; i++ )
                        {
                            T weight;
                            // declaring this const, helps the compiler to optimize (decision outside the for loop(s))
                            switch (method)
                            {
                                case(VPM::Gamma3):
                                    weight = interpolation_Kernel_Gamma3(new_pos[j] - old_pos[i], params.m_dx);
                                    break;
                                case(VPM::M4d):
                                    weight = interpolation_Kernel_M4d(new_pos[j] - old_pos[i], params.m_dx);
                                    break;
                                case(VPM::M5d):
                                    weight = interpolation_Kernel_M5d(new_pos[j] - old_pos[i], params.m_dx);
                                    break;
                                case(VPM::M6d):
                                    weight = interpolation_Kernel_M6d(new_pos[j] - old_pos[i], params.m_dx);
                                    break;
                            }
                            field[j] += weight * field_old[i];
                        }
                    }
                }
                break;
            //case(VPM::order6):
            //    std::cerr<<"    ---> remeshing order6"<<std::endl;
            //    for(unsigned int j=0; j<new_N; j++ )
            //    {
            //        for(unsigned int i=0; i<params.m_N; i++ )
            //        {
            //            T weight = interpolation_Kernel_order6(new_pos[j] - old_pos[i], params.m_dx);
            //            field[j] += weight * field_old[i];
            //        }
            //    }
            //    break;
            case(VPM::Gamma3_fast) : case(VPM::M4d_fast) : case(VPM::M5d_fast) : case(VPM::M6d_fast):
                {
                    int stenc = 2;
                    if (method == VPM::M6d_fast)
                    {
                        stenc = 5;
                    }
                    int num___ = 0;
                    for(unsigned int index_part=0; index_part<params.m_N; index_part++ )
                    {

                        if (g2p == RED2D_G2P_USE)
                        {
                            for (int i=0; i<m_grid_to_particles[index_part].size(); i++)
                            {
                                int particle_ind = m_grid_to_particles[index_part][i];
                                // declaring this const, helps the compiler to optimize (decision outside the for loop(s))
                                switch (method)
                                {
                                    case(VPM::Gamma3_fast):
                                        field[particle_ind] += interpolation_Kernel_Gamma3(new_pos[particle_ind] - old_pos[index_part], params.m_dx) * field_old[particle_ind];
                                        break;
                                    case(VPM::M4d_fast):
                                        field[particle_ind] += interpolation_Kernel_M4d(new_pos[particle_ind] - old_pos[index_part], params.m_dx) * field_old[particle_ind];
                                        break;
                                    case(VPM::M5d_fast):
                                        field[particle_ind] += interpolation_Kernel_M5d(new_pos[particle_ind] - old_pos[index_part], params.m_dx) * field_old[particle_ind];
                                        break;
                                    case(VPM::M6d_fast):
                                        field[particle_ind] += interpolation_Kernel_M6d(new_pos[particle_ind] - old_pos[index_part], params.m_dx) * field_old[particle_ind];
                                        break;
                                }
                            }
                        }
                        else
                        {
                            VPM::IPoint2d ind_min;
                            VPM::IPoint2d ind_max;
                            get_indices(params, stenc, old_pos[index_part], ind_min, ind_max);
                            for (unsigned int i=ind_min.x; i<=ind_max.x; i++)
                            {
                                for (unsigned int j=ind_min.y; j<=ind_max.y; j++)
                                {
                                    int index_reg = flat(i,j,int(params.m_num_px));
                                    // declaring this const, helps the compiler to optimize (decision outside the for loop(s))
                                    switch (method)
                                    {
                                        case(VPM::Gamma3_fast):
                                            field[index_reg] += interpolation_Kernel_Gamma3(new_pos[index_reg] - old_pos[index_part], params.m_dx) * field_old[index_part];
                                            break;
                                        case(VPM::M4d_fast):
                                            field[index_reg] += interpolation_Kernel_M4d(new_pos[index_reg] - old_pos[index_part], params.m_dx) * field_old[index_part];
                                            break;
                                        case(VPM::M5d_fast):
                                            field[index_reg] += interpolation_Kernel_M5d(new_pos[index_reg] - old_pos[index_part], params.m_dx) * field_old[index_part];
                                            break;
                                        case(VPM::M6d_fast):
                                            field[index_reg] += interpolation_Kernel_M6d(new_pos[index_reg] - old_pos[index_part], params.m_dx) * field_old[index_part];
                                            break;
                                    }
                                    if (g2p == RED2D_G2P_BUILD)
                                    {
                                        m_grid_to_particles[index_reg].push_back(index_part);
                                    }
                                }
                            }
                        }//if (g2p == RED2D_G2P_USE)
                    }
                }
                break;
            //case(VPM::order6_fast):
            //    {
            //        std::cerr<<"    ---> remeshing order6 fast"<<std::endl;
            //        int stenc = 5;
            //        for(unsigned int index_part=0; index_part<params.m_N; index_part++ )
            //        {

            //            VPM::IPoint2d ind_min;
            //            VPM::IPoint2d ind_max;
            //            get_indices(params, stenc, old_pos[index_part], ind_min, ind_max);

            //            for (unsigned int i=ind_min.x; i<=ind_max.x; i++)
            //            {
            //                for (unsigned int j=ind_min.y; j<=ind_max.y; j++)
            //                {
            //                    int index_reg = flat(i,j,int(params.m_num_px));
            //                    field[index_reg] += interpolation_Kernel_order6(new_pos[index_reg] - old_pos[index_part], params.m_dx) * field_old[index_part];
            //                }
            //            }
            //        }
            //    }
            //    break;
        }// end cases
        params.m_N = new_N;
    }

    template void Redistribute::redistribute<double>(
            const std::vector<Point2d> & old_pos,
            const std::vector<Point2d> & new_pos,
            Parameters & params,
            std::vector<double> & field,
            const REDISTRIBUTE2D g2p
            );

    template void Redistribute::redistribute<VPM::Point2d>(
            const std::vector<Point2d> & old_pos,
            const std::vector<Point2d> & new_pos,
            Parameters & params,
            std::vector<VPM::Point2d> & field,
            const REDISTRIBUTE2D g2p
            );

}

