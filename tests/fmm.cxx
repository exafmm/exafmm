#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>
#include "args.h"
#include "kernel.h"
#include "namespace.h"
#include "tree.h"
#include "verify.h"
#define CHECK_MONOPOLE

using namespace EXAFMM_NAMESPACE;

int main(int argc, char ** argv) {
  const real_t eps2 = 0.0;
  const complex_t wavek = complex_t(1.,.1) / real_t(2 * M_PI);
  Args args(argc, argv);

  Kernel kernel(args.P, eps2, wavek);
  logger::verbose = true;
  Verify verify;

  const int N = 50000;
  int level = 5;
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
  
  // key level offset (for global key)
  std::vector<int> keyLevelOffset(level+1, 0);
  for (l=0; l<=level; l++) keyLevelOffset[l] = ((1 << 3*l) - 1) / 7;

  // initialize M, L maps
  CoefMap M, L;
  
  // leaf level
  int key = -1;
  for (i=0; i<N; i++) {              // i: current particle index
    if (key != targetKeys[i]) {      // if previous key != current key (found a new nonempty cell)
      key = targetKeys[i];           // store the current key
      targetIndex2Key[level][targetNonEmpty[level]] = key;   // map nonempty counter -> new key
      targetKey2Index[level][key] = targetNonEmpty[level];   // map new key -> nonempty counter
      targetOffset[targetNonEmpty[level]] = i;               // record the first target index of this new cell
      targetNonEmpty[level]++;       // increase the counter of nonempty cells at leaf level
      L[key+keyLevelOffset[level]].resize(kernel.NTERM, 0.0);
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
      M[key+keyLevelOffset[level]].resize(kernel.NTERM, 0.0);
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
        L[key+keyLevelOffset[l-1]].resize(kernel.NTERM, 0.0);
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
        M[key+keyLevelOffset[l-1]].resize(kernel.NTERM, 0.0);
      }
    }
  }
  
  // level offset
  std::vector<int> targetLevelOffset(level+1, 0), sourceLevelOffset(level+1, 0);
  for (l=0; l<level; l++) {
    targetLevelOffset[l+1] = targetLevelOffset[l] + targetNonEmpty[l];
    sourceLevelOffset[l+1] = sourceLevelOffset[l] + sourceNonEmpty[l];
  }

#ifdef CHECK_LEVEL_OFFSET 
  for (l=0; l<level+1; l++) std::cout << targetNonEmpty[l] << " " 
                                      << targetLevelOffset[l] << std::endl;
  for (l=0; l<level+1; l++) std::cout << sourceNonEmpty[l] << " " 
                                      << sourceLevelOffset[l] << std::endl;
  for (l=0; l<level+1; l++) std::cout << keyLevelOffset[l] << " " << std::endl;
#endif
  
  int ic, jc, j;
  int ni, nj;
  vec3 Xi, Xj, XI, XJ;

  // P2M
  for (jc=0; jc<sourceNonEmpty[level]; jc++) {              // jc: nonempty source cell loop counter
    key = sourceIndex2Key[level][jc];
    getX(Xj, key, level);
    Source* Sj = &sources[sourceOffset[jc]];
    nj = sourceOffset[jc+1] - sourceOffset[jc];
    kernel.P2M(Xj, M[key+keyLevelOffset[level]], Sj, nj);
  }

  // M2M
  for (l=level; l>0; l--) {                   // l: current level
    for (jc=0; jc<sourceNonEmpty[l]; jc++) {  // jc: nonempty source cell index in level l
      key = sourceIndex2Key[l][jc];           // key: jc's Morton key
      getX(Xj, key, l);                       // Xj: child cell center
      getX(XJ, key/8, l-1);                   // XJ: parent cell center
      kernel.M2M(XJ, M[key/8+keyLevelOffset[l-1]], Xj, M[key+keyLevelOffset[l]]);
    }
  }
  
  // check monopole of root
#ifdef CHECK_MONOPOLE
  std::cout << "monopole of root: " << M[0][0] << std::endl;
#endif

  // M2L
  for (l=2; l<=level; l++) {                // l: current level
    for (ic=0; ic<targetNonEmpty[l]; ic++) {   // i: nonempty target cell index in level l
      getIX(iX, targetIndex2Key[l][ic]);     // i -> key -> iX
      for (jc=0; jc<sourceNonEmpty[l]; jc++) { // j: nonempty source cell index in level l
        getIX(jX, sourceIndex2Key[l][jc]);   // j -> key -> jX
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
  
  /*

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
  */
}
