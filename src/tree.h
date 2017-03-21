#ifndef tree_h
#define tree_h
#include <cmath>
#include "namespace.h"
#include "types.h"

namespace EXAFMM_NAMESPACE {
  uint64_t getKey(ivec3 iX, int level);
  void getIX(ivec3 & iX, uint64_t key);
  void getX(vec3 & X, uint64_t key, int level);
  void sortKey(std::vector<int> keys, std::vector<int> & permutation);
  
  template<typename T>
  void permutate(T & vec, std::vector<int> permutation) {
    T vec2(vec.size());
    for (int i=0; i<vec.size(); i++) vec2[permutation[i]] = vec[i];
    for (int i=0; i<vec.size(); i++) vec[i] = vec2[i];
  }

  inline bool isNeighbor(ivec3 iX, ivec3 jX) {
    return (std::abs(iX[0]-jX[0]) <= 1) && 
           (std::abs(iX[1]-jX[1]) <= 1) && 
           (std::abs(iX[2]-jX[2]) <= 1); 
  }

}
#endif
