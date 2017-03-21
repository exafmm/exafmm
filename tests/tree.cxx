#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>
#include "args.h"
#include "namespace.h"
#include "tree.h"

#define CHECK_MONOPOLE

using namespace EXAFMM_NAMESPACE;

int main(int argc, char ** argv) {
  const int N = 10000;
  int level = 4;
  int nx = 1 << level;          // number of leafs in each dimension
  int nleaf = nx * nx * nx;     // number of leafs if full
  int i, l, d;                  // iterator in loops

  Targets targets(N);
  Sources sources(N);
  std::vector<int> targetKeys(N), sourceKeys(N);
  ivec3 iX, jX;

  // generate rand X and key for targets and sources
  for (i=0; i<N; i++) {
    for (d=0; d<3; d++) {
      targets[i].X[d] = drand48();
      sources[i].X[d] = drand48();
      iX[d] = (int) (targets[i].X[d] * nx);
      jX[d] = (int) (sources[i].X[d] * nx);
    }
    targets[i].F = 0;
    sources[i].Q = 1.0/N;
    targetKeys[i] = getKey(iX, level);
    sourceKeys[i] = getKey(jX, level);
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
  
  std::vector<real_t> M(sourceNonEmpty[level]+sourceLevelOffset[level], 0.0);
  std::vector<real_t> L(targetNonEmpty[level]+targetLevelOffset[level], 0.0);

  // P2M
  int j;
   for (i=0; i<sourceNonEmpty[level]; i++) {              // i: nonempty source cell index
    for (j=sourceOffset[i]; j<sourceOffset[i+1]; j++) {   // j: source index
      M[i+sourceLevelOffset[level]] += sources[j].Q;
    }
  }

  // M2M
  for (l=level; l>0; l--) {                 // l: current level
    for (i=0; i<sourceNonEmpty[l]; i++) {   // i: nonempty source cell index in level l
      key = sourceIndex2Key[l][i];          // key: i's Morton key
      M[sourceKey2Index[l-1][key/8]+sourceLevelOffset[l-1]] += M[i+sourceLevelOffset[l]];
    }
  }
  
  // check monopole of root
#ifdef CHECK_MONOPOLE
  std::cout << "monopole of root: " << M[0] << std::endl;
  return 0;
#endif

  // M2L
  int count=0;
  for (l=2; l<=level; l++) {                // l: current level
    nx = 1 << l;                            // nx: number of cells in one-axis in level l
    for (i=0; i<targetNonEmpty[l]; i++) {   // i: nonempty target cell index in level l
      getIX(iX, targetIndex2Key[l][i]);     // i -> key -> iX
      for (j=0; j<sourceNonEmpty[l]; j++) { // j: nonempty source cell index in level l
        getIX(jX, sourceIndex2Key[l][j]);   // j -> key -> jX
        if (isNeighbor(iX/2, jX/2)) {       // if i,j's parents are neighbor
          if (!isNeighbor(iX, jX)) {        // and i,j are non-neighbor
            count++;
            real_t dx = (real_t) (iX[0] - jX[0]) / nx;
            real_t dy = (real_t) (iX[1] - jX[1]) / nx;
            real_t dz = (real_t) (iX[2] - jX[2]) / nx;
            real_t r = std::sqrt(dx*dx + dy*dy + dz*dz);
            L[i+targetLevelOffset[l]] += M[j+sourceLevelOffset[l]] / r;
          }         
        }
      }
    }
  }

  // L2L
  for (l=3; l<=level; l++) {                // l: current level
    for (i=0; i<targetNonEmpty[l]; i++) {   // i: nonempty target cell index in level l
      key = targetIndex2Key[l][i];          // key: i's Morton key 
      L[i+targetLevelOffset[l]] += L[targetKey2Index[l-1][key/8]+targetLevelOffset[l-1]];
    }
  }

  // L2P
  for (i=0; i<targetNonEmpty[level]; i++) {
    for (j=targetOffset[i]; j<targetOffset[i+1]; j++) {
      targets[j].F += L[i+targetLevelOffset[level]];
    }
  }

  // P2P
  for (int ic=0; ic<targetNonEmpty[level]; ic++) {
    getIX(iX, targetIndex2Key[level][ic]);
    for (int jc=0; jc<sourceNonEmpty[level]; jc++) {
      getIX(jX, sourceIndex2Key[level][jc]);
      if (isNeighbor(iX, jX)) {
        for (i=targetOffset[ic]; i<targetOffset[ic+1]; i++) {
          for (j=sourceOffset[jc]; j<sourceOffset[jc+1]; j++) {
          vec3 dX = targets[i].X - sources[j].X;
          real_t r = std::sqrt(norm(dX));
          if (r!=0) targets[i].F += sources[j].Q / r;
          }
        }
      }
    }
  }

#ifdef PRINT_RESULT
  // check answer
  for (i=0; i<N; i++) {
    real_t Fi = 0;
    for (j=0; j<N; j++) {
      vec3 dX = targets[i].X - sources[j].X;
      real_t r = std::sqrt(norm(dX));
      if (r!=0) Fi += sources[j].Q / r;
    }
    std::cout << i << " " << targets[i].F << " " << Fi << std::endl;
  }
#endif
}
