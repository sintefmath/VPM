#pragma once
#include "Point3d.hpp"
#include "FlattenArrayIndex.hpp"
#include <vector>

/**
 * (u,v)
 * @returns dv/dx - du/dy
 */
inline VPM::Point3d operator_curl(
        const double fx_i_j_km, const double fx_i_j_kp,
        const double fx_i_jm_k, const double fx_i_jp_k,
        const double fy_im_j_k, const double fy_ip_j_k,
        const double fy_i_j_km, const double fy_i_j_kp,
        const double fz_im_j_k, const double fz_ip_j_k,
        const double fz_i_jm_k, const double fz_i_jp_k,
        const double dx, const double dy, const double dz
        )
{
    return VPM::Point3d(
            // +dpz/dy - dpy/dz
            (fz_i_jp_k - fz_i_jm_k)/(2.*dy) - (fy_i_j_kp - fy_i_j_km)/(2.*dz),
            // +dpx/dz - dpz/dx
            (fx_i_j_kp - fx_i_j_km)/(2.*dz) - (fz_ip_j_k - fz_im_j_k)/(2.*dx),
            // +dpy/dx - dpx/dy
            (fy_ip_j_k - fy_im_j_k)/(2.*dx) - (fx_i_jp_k - fx_i_jm_k)/(2.*dy)
            );
}


inline std::vector<VPM::Point3d> calculate_curl(
        const std::vector<VPM::Point3d> & velocity,
        const unsigned int num_px,
        const unsigned int num_py,
        const unsigned int num_pz,
        const double dx,
        const double dy,
        const double dz
        )
{
    std::vector<VPM::Point3d> curl;
    curl.resize(velocity.size(), VPM::Point3d(0));

    for (int k=0; k<num_pz; k++)
    {
        int km = std::max(k-1,0);
        int kp = std::min(k+1,int(num_pz-1));
        for (int j=0; j<num_py; j++)
        {
            int jm = std::max(j-1,0);
            int jp = std::min(j+1,int(num_py-1));
            for (int i=0; i<num_px; i++)
            {
                int im = std::max(i-1,0);
                int ip = std::min(i+1,int(num_px-1));

                int i_j_k = flat(i,j,k, int(num_px), int(num_py));

                int im_j_k = flat(im,j,k, int(num_px), int(num_py));
                int ip_j_k = flat(ip,j,k, int(num_px), int(num_py));
                int i_jm_k = flat(i,jm,k, int(num_px), int(num_py));
                int i_jp_k = flat(i,jp,k, int(num_px), int(num_py));
                int i_j_km = flat(i,j,km, int(num_px), int(num_py));
                int i_j_kp = flat(i,j,kp, int(num_px), int(num_py));

                curl[i_j_k] = operator_curl(
                    velocity[i_j_km].x, velocity[i_j_kp].x,
                    velocity[i_jm_k].x, velocity[i_jp_k].x,
                    velocity[im_j_k].y, velocity[ip_j_k].y,
                    velocity[i_j_km].y, velocity[i_j_kp].y,
                    velocity[im_j_k].z, velocity[ip_j_k].z,
                    velocity[i_jm_k].z, velocity[i_jp_k].z,
                    dx, dy, dz
                    );

            }
        }
    }
    return curl;
}

