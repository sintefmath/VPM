#include "Split_Source.hpp"
#include "Utilities.hpp"

namespace VPM
{

    Split_Source::Split_Source(
            )
    {
    };
    Split_Source::~Split_Source()
    {
    };

//    void Split_Source::setVelocityFields(
//            const std::vector<Point3d> & velocity,
//            const std::vector<Point3d> & background_velocity,
//            const std::shared_ptr<Parameters> params
//            )
//    {
//        //m_velocity = velocity;
//        //m_background_velocity = background_velocity;
//        m_total_velocity.resize(velocity.size());
//        for(unsigned int i=0; i<velocity.size(); i++ )
//        {
//            m_total_velocity[i] = velocity[i] + background_velocity[i];
//        }
//        m_curl_background_velocity = calculate_curl(background_velocity, params->m_num_px, params->m_num_py, params->m_dx, params->m_dy);
//    }
//
//
//
//    std::vector<double> Split_Source::step(
//            const std::vector<double> & omega,
//            const std::shared_ptr<Parameters> params
//            )
//    {
//
//        // return -(u+v)\nabla curl(v) + \nu laplace curl(v)
//
//        std::vector<double> ret_source;
//        ret_source.resize(omega.size());
//
//        double nu_lambda = params->m_nu/(params->m_dx*params->m_dx);
//
//        for (int j=0; j<params->m_num_py; j++)
//        {
//            for (int i=0; i<params->m_num_px; i++)
//            {
//                int im = std::max(i-1,0);
//                int ip = std::min(i+1,int(params->m_num_px-1));
//                int jm = std::max(j-1,0);
//                int jp = std::min(j+1,int(params->m_num_py-1));
//                int i_j   = i +int(params->m_num_px)*j;
//                int i_jm  = i +int(params->m_num_px)*jm;
//                int i_jp  = i +int(params->m_num_px)*jp;
//                int im_j  = im+int(params->m_num_px)*j;
//                //int im_jm = im+int(params->m_num_px)*jm;
//                //int im_jp = im+int(params->m_num_px)*jp;
//                int ip_j  = ip+int(params->m_num_px)*j;
//                //int ip_jm = ip+int(params->m_num_px)*jm;
//                //int ip_jp = ip+int(params->m_num_px)*jp;
//                ret_source[i_j] =
//                    - dot(
//                        m_total_velocity[i_j],
//                        Point3d(
//                            // central difference for curl(v)_x
//                            (m_curl_background_velocity[ip_j] - m_curl_background_velocity[im_j])/params->m_dx,
//                            // central difference for curl(v)_y
//                            (m_curl_background_velocity[i_jp] - m_curl_background_velocity[i_jm])/params->m_dy 
//                            )
//                        )
//                    +
//                    nu_lambda*
//                    // cur(v)_xx
//                    ((m_curl_background_velocity[ip_j] - 2*m_curl_background_velocity[i_j] + m_curl_background_velocity[im_j])
//                    +
//                    // cur(v)_yy
//                     (m_curl_background_velocity[i_jp] - 2*m_curl_background_velocity[i_j] + m_curl_background_velocity[i_jm]))
//                    ;
//            }
//        }
//        return ret_source;
//    }

}
