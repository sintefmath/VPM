#pragma once
#include "Point2d.hpp"
#include <vector>

/**
 * (u,v)
 * @returns dv/dx - du/dy
 */
inline double operator_curl(const double Ux_i_jp, const double Ux_i_jm, const double Uy_ip_j, const double Uy_im_j, const double dx, const double dy)
{
    return (Uy_ip_j - Uy_im_j)/(2.*dx) - (Ux_i_jp - Ux_i_jm)/(2.*dy);
}


inline std::vector<double> calculate_curl(
        const std::vector<VPM::Point2d> & velocity,
        const unsigned int num_px,
        const unsigned int num_py,
        const double dx,
        const double dy
        )
{
    std::vector<double> curl;
    curl.resize(velocity.size(), double(0.));

    for (int j=0; j<num_py; j++)
    {
        for (int i=0; i<num_px; i++)
        {
            int im = std::max(i-1,0);
            int ip = std::min(i+1,int(num_px-1));
            int jm = std::max(j-1,0);
            int jp = std::min(j+1,int(num_py-1));
            int i_j   = i +int(num_px)*j;
            int i_jm  = i +int(num_px)*jm;
            int i_jp  = i +int(num_px)*jp;
            int im_j  = im+int(num_px)*j;
            int ip_j  = ip+int(num_px)*j;
            curl[i_j] = operator_curl(velocity[i_jp].x, velocity[i_jm].x, velocity[ip_j].y, velocity[im_j].y, dx, dy);
        }
    }
    return curl;
}

