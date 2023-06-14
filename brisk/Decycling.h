#ifndef _DECYCLING
#define _DECYCLING

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <cstdint>



using namespace std;



class DecyclingSet {
  private:
    const uint k;
    const double unit;
    vector<double> coef;
    const double eps = 0.000001;
    double computeR(uint64_t seq);

  public:
    DecyclingSet(uint k);
    ~DecyclingSet() {}
    bool mem(uint64_t seq);
    uint memDouble(uint64_t seq);
};

#endif
