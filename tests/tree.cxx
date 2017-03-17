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
  int i, l, d;                  // iterator in loops

  Targets targets(N);
  Sources sources(N);
  std::vector<int> targetKeys(N), sourceKeys(N);
  ivec3 iX, iX2;

  // generate rand X and key for targets and sources
  for (i=0; i<N; i++) {
    for (d=0; d<3; d++) {
      targets[i].X[d] = drand48();
      sources[i].X[d] = drand48();
      iX[d] = (int) (targets[i].X[d] * nx);
      iX2[d] = (int) (sources[i].X[d] * nx);
    }
    targets[i].F = 0;
    sources[i].Q = 1/N;
    targetKeys[i] = getKey(iX, level);
    sourceKeys[i] = getKey(iX2, level);
  }

  // sort targets and sources based on Morton keys
  std::vector<int> targetPermutation(N), sourcePermutation(N);
  sortKey(targetKeys, targetPermutation);
  sortKey(sourceKeys, sourcePermutation);
  
  // permutate targets, sources and keys
  permutate<Targets>(targets, targetPermutation);  
  permutate<Sources>(sources, sourcePermutation);  
  permutate< std::vector<int> >(targetKeys, targetPermutation);
  permutate< std::vector<int> >(sourceKeys, sourcePermutation);
  
#if 0
  for (i=0; i<N; i++) {
    std::cout << i << " " << targetKeys[i] << " " << targets[i].X << std::endl;
    std::cout << i << " " << sourceKeys[i] << " " << sources[i].X << std::endl;
  }
#endif 

  // calculate offset vectors for adaptive tree
  std::vector<int> targetNonEmpty(level+1, 0), sourceNonEmpty(level+1, 0);  // num of non-empty cells at each level
  std::vector<int> targetOffset(nleaf), sourceOffset(nleaf);  // leaf morton key -> particle index offset
  std::vector< std::vector<int> > targetIndex2Key( level+1,    // index2Key[level][index] -> Morton key
                                                   std::vector<int>(nleaf, 0) );
  std::vector< std::vector<int> > sourceIndex2Key( level+1,
                                                   std::vector<int>(nleaf, 0) );
  std::vector< std::vector<int> > targetKey2Index( level+1,    // key2Index[level][key] -> index at non-empty cells' list
                                                   std::vector<int>(nleaf, 0) );
  std::vector< std::vector<int> > sourceKey2Index( level+1,
                                                   std::vector<int>(nleaf, 0) );
  
  // leaf level
  int key = -1;
  for (i=0; i<N; i++) {              // i: current particle index
    if (key != targetKeys[i]) {      // if previous key != current key (found a new nonempty cell)
      key = targetKeys[i];           // store the current key
      targetIndex2Key[level][targetNonEmpty[level]] = key;   // map nonempty counter -> new key
      targetKey2Index[level][key] = targetNonEmpty[level];   // map new key -> nonempty counter
      targetOffset[targetNonEmpty[level]] = i;               // record the first target index of this new cell
      targetNonEmpty[level]++;       // increase the counter of nonempty cells at leaf level
    }
  }
  targetOffset[targetNonEmpty[level]] = N;                   // for the last leaf cell

  key = -1;
  for (i=0; i<N; i++) {              // i: current particle index
    if (key != sourceKeys[i]) {      // if previous key != current key (found a new nonempty cell)
      key = sourceKeys[i];           // store the current key
      sourceIndex2Key[level][sourceNonEmpty[level]] = key;   // map nonempty counter -> new key
      sourceKey2Index[level][key] = sourceNonEmpty[level];   // map new key -> nonempty counter
      sourceOffset[sourceNonEmpty[level]] = i;               // record the first target index of this new cell
      sourceNonEmpty[level]++;       // increase the counter of nonempty cells at leaf level
    }
  }
  sourceOffset[sourceNonEmpty[level]] = N;                   // for the last leaf cell
  
  // non-leaf level
  for (l=level; l>0; l--) {          // l: current level, l-1: parent level
    key = -1;                        // reset key
    for (i=0; i<targetNonEmpty[l]; i++) {                  // i: nonempty cell at level l
      if (key != targetIndex2Key[l][i]/8) {                // if previous key != i's parent key 
        key = targetIndex2Key[l][i]/8;                     // store i's parent key
        targetIndex2Key[l-1][targetNonEmpty[l-1]] = key;   // map nonempty counter -> new parent key 
        targetKey2Index[l-1][key] = targetNonEmpty[l-1];   // map new parent key -> nonempty counter
        targetNonEmpty[l-1]++;                             // 
      }
    }
  }
  for (l=level; l>0; l--) {          // l: current level, l-1: parent level
    key = -1;                        // reset key
    for (i=0; i<sourceNonEmpty[l]; i++) {                  // i: nonempty cell at level l
      if (key != sourceIndex2Key[l][i]/8) {                // if previous key != i's parent key 
        key = sourceIndex2Key[l][i]/8;                     // store i's parent key
        sourceIndex2Key[l-1][sourceNonEmpty[l-1]] = key;   // map nonempty counter -> new parent key 
        sourceKey2Index[l-1][key] = sourceNonEmpty[l-1];   // map new parent key -> nonempty counter
        sourceNonEmpty[l-1]++;                             // 
      }
    }
  }

  // level offset
  std::vector<int> targetLevelOffset(level+1, 0), sourceLevelOffset(level+1, 0);
  for (l=0; l<level; l++) {
    targetLevelOffset[l+1] = targetLevelOffset[l] + targetNonEmpty[l];
    sourceLevelOffset[l+1] = sourceLevelOffset[l] + sourceNonEmpty[l];
  }

#if 0
  for (l=0; l<level+1; l++) std::cout << targetNonEmpty[l] << " " 
                                      << targetLevelOffset[l] << std::endl;
  for (l=0; l<level+1; l++) std::cout << sourceNonEmpty[l] << " " 
                                      << sourceLevelOffset[l] << std::endl;
#endif

}

