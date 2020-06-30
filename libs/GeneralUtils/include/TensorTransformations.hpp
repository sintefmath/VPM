#pragma once
#include <FlattenArrayIndex.hpp>
#include <vector>
#include <cmath>
#include <assert.h>

inline void cpu_crop_permute_convert_subtract(std::vector<std::vector<unsigned char>> input, std::vector<int> insize_wh, std::vector<float> & output, std::vector<int> pad_wh, std::vector<int> outsize_wh, std::vector<float> mean_rgb, int j, int k)
{
    const std::vector<unsigned char> px = input[ flat(pad_wh[0]+j,pad_wh[1]+k, insize_wh[0]) ];

    output[flat(j,k,0, outsize_wh[0],outsize_wh[1])] = float(px[0])-mean_rgb[0];
    output[flat(j,k,1, outsize_wh[0],outsize_wh[1])] = float(px[1])-mean_rgb[1];
    output[flat(j,k,2, outsize_wh[0],outsize_wh[1])] = float(px[2])-mean_rgb[2];
}

template<typename type>
inline void permuteIK(std::vector<type> & tensor, int & I, int & J, int & K)
{

    assert(I*J*K==tensor.size());

    std::vector<type> tmp = tensor;

    for (int k = 0; k < K; k++)
    {
        for (int j = 0; j < J; j++)
        {
            for (int i = 0; i < I; i++)
            {
                tensor[flat(j,k,i, J,K)] = tmp[flat(i,j,k, I, J)];
            }
        }
    }

    int tmp_I=I;
    I = K;
    K = tmp_I;

}

template<typename type>
inline void flipI(std::vector<type> & tensor, int I, int J, int K)
{

    assert(I*J*K==tensor.size());

    std::vector<type> tmp = tensor;

    for (int i = 0; i < I; i++)
    {
        for (int j = 0; j < J; j++)
        {
            for (int k = 0; k < K; k++)
            {
                tensor[flat(I-1-i,j,k, I,J)] = tmp[flat(i,j,k, I,J)];
            }
        }
    }

}

template<typename type>
inline void subtract(std::vector<type> & tensor, type val)
{

    for (int i = 0; i < tensor.size(); i++)
    {
        tensor[i] -= val;
    }

}

template<typename type, typename type2>
inline void convert(const std::vector<type> & tensor, std::vector<type2> & conv_tensor)
{

    conv_tensor.resize(tensor.size());
    for (int i = 0; i < tensor.size(); i++)
    {
        conv_tensor[i] = static_cast<type2>(tensor[i]);
    }

}

template<typename type>
inline void flipJ(std::vector<type> & tensor, int I, int J, int K)
{

    assert(I*J*K==tensor.size());

    std::vector<type> tmp = tensor;

    for (int i = 0; i < I; i++)
    {
        for (int j = 0; j < J; j++)
        {
            for (int k = 0; k < K; k++)
            {
                tensor[flat(i,J-1-j,k, I,J)] = tmp[flat(i,j,k, I,J)];
            }
        }
    }

}

template<typename type>
inline void flipK(std::vector<type> & tensor, int I, int J, int K)
{

    assert(I*J*K==tensor.size());

    std::vector<type> tmp = tensor;

    for (int i = 0; i < I; i++)
    {
        for (int j = 0; j < J; j++)
        {
            for (int k = 0; k < K; k++)
            {
                tensor[flat(i,j,K-1-k, I,J)] = tmp[flat(i,j,k, I,J)];
            }
        }
    }

}

template<typename type>
inline void centercropIJ(std::vector<type> & tensor, int I_old, int J_old, int K, int I_new, int J_new)
{

    assert(I_old*J_old*K==tensor.size());
    assert(I_new<I_old);
    assert(J_new<J_old);

    int pad_I = std::floor(float(I_old - I_new)/2.f);
    int pad_J = std::floor(float(J_old - J_new)/2.f);

    std::vector<type> tmp = tensor;
    tensor.resize(I_new*J_new*K);

    for (int i = 0; i < I_new; i++)
    {
        for (int j = 0; j < J_new; j++)
        {
            for (int k = 0; k < K; k++)
            {
                tensor[flat(i,j,k, I_new,J_new)] = tmp[flat(pad_I+i,pad_J+j,k, I_old,J_old)];
            }
        }
    }
}

template<typename type>
inline void centercropJK(std::vector<type> & tensor, int I, int J_old, int K_old, int J_new, int K_new)
{

    assert(I*J_old*K_old==tensor.size());
    assert(J_new<J_old);
    assert(K_new<K_old);

    int pad_J = std::floor(float(J_old - J_new)/2.f);
    int pad_K = std::floor(float(K_old - K_new)/2.f);

    std::vector<type> tmp = tensor;
    tensor.resize(I*J_new*K_new);

    for (int i = 0; i < I; i++)
    {
        for (int j = 0; j < J_new; j++)
        {
            for (int k = 0; k < K_new; k++)
            {
                tensor[flat(i,j,k, I,J_new)] = tmp[flat(i,pad_J+j,pad_K+k, I,J_old)];
            }
        }
    }
}
