#pragma once

inline float linear_ab_to_cd(int x, int a, int b, int c, int d)
{
    return ( ((float)d - (float)c) / ((float)b - (float)a) * ((float)x - (float)a) + (float)c );
}

inline float float_linear_ab_to_cd(float x, float a, float b, float c, float d)
{
    return ( (d - c) / (b - a) * (x - a) + c );
}

// Round a / b to nearest higher integer value
inline int iDivUp(int a, int b)
{
   return (a % b != 0) ? (a / b + 1) : (a / b);
}
