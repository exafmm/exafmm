#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <vector>
#include "args.h"
#include "kernel.h"
#include "namespace.h"
#include "verify.h"

using namespace EXAFMM_NAMESPACE;

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

// permutate a vector according to its permutation
template<typename T>
void permutate(T & vec, std::vector<int> permutation) {
  T vec2(vec.size());
  for (int i=0; i<vec.size(); i++) vec2[permutation[i]] = vec[i];
  for (int i=0; i<vec.size(); i++) vec[i] = vec2[i];
}


int main(int argc, char ** argv) {
  const int N = 100;
  int level = 3;
  int nx = 1 << level;          // number of leafs in each dimension
  int nleaf = nx * nx * nx;     // number of leafs if full

  Targets targets(N), targets2(N);
  Sources sources(N);
  std::vector<int> targetCells(N), sourceCells(N);
  ivec3 iX, iX2;

  // generate rand X and key for targets and sources
  for (int i=0; i<N; i++) {
    for (int d=0; d<3; d++) {
      targets[i].X[d] = drand48();
      sources[i].X[d] = drand48();
      iX[d] = (int) (targets[i].X[d] * nx);
      iX2[d] = (int) (sources[i].X[d] * nx);
    }
    targets[i].F = 0;
    sources[i].Q = 1/N;
    targetCells[i] = getKey(iX, level);
    sourceCells[i] = getKey(iX2, level);
  }

  // sort targets and sources based on Morton keys
  std::vector<int> targetPermutation(N), sourcePermutation(N);
  sortKey(targetCells, targetPermutation);
  sortKey(sourceCells, sourcePermutation);
  
  // permutate targets, sources and keys
  permutate<Targets>(targets, targetPermutation);  
  permutate<Sources>(sources, sourcePermutation);  
  permutate< std::vector<int> >(targetCells, targetPermutation);
  permutate< std::vector<int> >(sourceCells, sourcePermutation);
  
#if 0
  for (int i=0; i<N; i++) {
    std::cout << i << " " << targetCells[i] << " " << targets[i].X << std::endl;
    std::cout << i << " " << sourceCells[i] << " " << sources[i].X << std::endl;
  }
#endif 

  // calculate offset vectors for adaptive tree
  std::vector<int> ncells(level+1);  // num of non-empty cells at each level
  std::vector<int> targetOffset(nleaf), sourceOffset(nleaf);  // leaf morton key --> particle index offset
   
}

