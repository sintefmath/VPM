#pragma once

template<typename T>
inline T clamp(T val, T val_min, T val_max)
{
    return std::max(val_min, std::min(val, val_max));
};
