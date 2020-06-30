#pragma once
#include <algorithm>

template<typename type, typename type2>//typically unsigned int or int...
#ifdef __CUDACC__
__host__ __device__
#endif
inline type flat(type i, type j, type k, type l, type m, type2 I, type2 J, type2 K, type2 L)
{
    return (i + I*(j + J*(k + K*(l + L*m))));
}

template<typename type, typename type2>//typically unsigned int or int...
#ifdef __CUDACC__
__host__ __device__
#endif
inline type flat(type i, type j, type k, type l, type2 I, type2 J, type2 K)
{
    return (i + I*(j + J*(k + K*l)));
}

template<typename type, typename type2>//typically unsigned int or int...
#ifdef __CUDACC__
__host__ __device__
#endif
inline type flat(type i, type j, type k, type2 I, type2 J)
{
    return (i + I*(j + J*k));
}

template<typename type, typename type2>//typically unsigned int or int...
#ifdef __CUDACC__
__host__ __device__
#endif
inline type flat(type i, type j, type2 I)
{
    return (i + I*j);
}

template<typename type>//typically unsigned int or int...
#ifdef __CUDACC__
__host__ __device__
#endif
inline type flat(type i, type j, type k, type l, type m, type I, type J, type K, type L)
{
    return (i + I*(j + J*(k + K*(l + L*m))));
}

template<typename type>//typically unsigned int or int...
#ifdef __CUDACC__
__host__ __device__
#endif
inline type flat(type i, type j, type k, type l, type I, type J, type K)
{
    return (i + I*(j + J*(k + K*l)));
}

template<typename type>//typically unsigned int or int...
#ifdef __CUDACC__
__host__ __device__
#endif
inline type flat(type i, type j, type k, type I, type J)
{
    return (i + I*(j + J*k));
}

template<typename type>//typically unsigned int or int...
#ifdef __CUDACC__
__host__ __device__
#endif
inline type flat(type i, type j, type I)
{
    return (i + I*j);
}

template<typename type>//typically unsigned int or int...
#ifdef __CUDACC__
__host__ __device__
#endif
inline type flatZeroNeumann(type i, type j, type k, type l, type m, type I, type J, type K, type L, type M)
{
    type a = std::min(std::max(i,0),I-1);
    type b = std::min(std::max(j,0),J-1);
    type c = std::min(std::max(k,0),K-1);
    type d = std::min(std::max(l,0),L-1);
    type e = std::min(std::max(m,0),M-1);
    return flat(a,b,c,d,e, I,J,K,L);
}

template<typename type>//typically unsigned int or int...
#ifdef __CUDACC__
__host__ __device__
#endif
inline type flatZeroNeumann(type i, type j, type k, type l, type I, type J, type K, type L)
{
    type a = std::min(std::max(i,0),I-1);
    type b = std::min(std::max(j,0),J-1);
    type c = std::min(std::max(k,0),K-1);
    type d = std::min(std::max(l,0),L-1);
    return flat(a,b,c,d, I,J,K);
}

template<typename type>//typically unsigned int or int...
#ifdef __CUDACC__
__host__ __device__
#endif
inline type flatZeroNeumann(type i, type j, type k, type I, type J, type K)
{
    type a = std::min(std::max(i,0),I-1);
    type b = std::min(std::max(j,0),J-1);
    type c = std::min(std::max(k,0),K-1);
    return flat(a,b,c, I,J);
}

template<typename type>//typically unsigned int or int...
#ifdef __CUDACC__
__host__ __device__
#endif
inline type flatZeroNeumann(type i, type j, type I, type J)
{
    type a = std::min(std::max(i,0),I-1);
    type b = std::min(std::max(j,0),J-1);
    return flat(a,b, I);
}

