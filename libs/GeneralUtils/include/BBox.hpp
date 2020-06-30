#pragma once

template <typename vecX>
struct BBox
{
    BBox( const vecX& min, const vecX& max )
        : m_min( min ), m_max( max )
    {}
    vecX m_min;
    vecX m_max;
};

template <>
struct BBox<float>
{
    BBox( const float& min, const float& max )
        : m_min( min ), m_max( max )
    {
    }
    float m_min;
    float m_max;
};

