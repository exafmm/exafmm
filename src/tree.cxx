#include <algorithm>
#include "tree.h"

namespace EXAFMM_NAMESPACE {
  // interleave iX into Morton key
  uint64_t getKey(ivec3 iX, int level) {
    // uint64_t index = ((1 << 3 * level) - 1) / 7;          // Level offset
    uint64_t key = 0;
    for (int l=0; l<level; l++) {
      for (int d=0; d<3; d++) { 
        key += (iX[d] & 1) << (3 * l + d);
        iX[d] >>= 1;
      }
    } 
    return key;
  }

  // deinterleave Morton key into iX
  void getIX(ivec3 & iX, uint64_t key) {
    iX = 0;
    int l = 0;
    while (key > 0) {
      for (int d=0; d<3; d++) {
        iX[d] += (key & 1) << l;
        key >>= 1;
      }
      l++;
    }
  }

  // compute cell's center given key and level
  void getX(vec3 & X, uint64_t key, int level) {
    ivec3 iX;
    getIX(iX, key);
    int nx = 1 << level;
    for (int d=0; d<3; d++) {
      X[d] = 1.0 / nx * (iX[d] + 0.5);
    }
  }

  // bucket sort, return permutation indices
  void sortKey(std::vector<int> keys, std::vector<int> & permutation) {
    int iMax = * std::max_element(keys.begin(), keys.end());
    std::vector<int> buckets(iMax+1, 0);
    int i, inew;
    for (i=0; i<keys.size(); i++) buckets[keys[i]]++;
    for (i=1; i<buckets.size(); i++) buckets[i] += buckets[i-1];
    for (i=keys.size()-1; i>=0; i--) {
      buckets[keys[i]]--;
      inew = buckets[keys[i]];
      permutation[i] = inew;
    }
  }
}
