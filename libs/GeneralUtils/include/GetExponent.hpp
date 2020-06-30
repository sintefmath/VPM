#include <cmath>

inline int getexp10(float arg)
{
   return std::floor(std::log10(std::fabs(arg) ));
}
