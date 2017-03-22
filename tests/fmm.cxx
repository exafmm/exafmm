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

  const int N = 1000;
  int level = 3;
  int nx = 1 << level;          // number of leafs in each dimension
  int nleaf = nx * nx * nx;     // number of leafs if full
  int i, l, d;                  // iterator in loops

  Targets targets(N), targets2(N);
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
  
  targets2 = targets;

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
  
  // key level offset (for global key)
  std::vector<int> keyLevelOffset(level+1, 0);
  for (l=0; l<=level; l++) keyLevelOffset[l] = ((1 << 3*l) - 1) / 7;

  // initialize M, L maps
  CoefMap M, L;
  
  // leaf level
  uint64_t key = -1;
  for (i=0; i<N; i++) {              // i: current particle index
    if (key != targetKeys[i]) {      // if previous key != current key (found a new nonempty cell)
      key = targetKeys[i];           // store the current key
      targetIndex2Key[level][targetNonEmpty[level]] = key;   // map nonempty counter -> new key
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

  int ic, jc, j;
  int ni, nj;
  vec3 Xi, Xj, XI, XJ;
  Target* Ti;
  Source* Sj;

  // P2M
  for (jc=0; jc<sourceNonEmpty[level]; jc++) {              // jc: nonempty source cell loop counter
    key = sourceIndex2Key[level][jc];
    getX(Xj, key, level);
    Sj = &sources[sourceOffset[jc]];
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
  uint64_t keyi, keyj;
  for (l=2; l<=level; l++) {                  // l: current level
    for (ic=0; ic<targetNonEmpty[l]; ic++) {  // ic: nonempty target cell index in level l
      keyi = targetIndex2Key[l][ic];
      getIX(iX, keyi);                        // ic -> keyi -> iX
      for (jc=0; jc<sourceNonEmpty[l]; jc++) {// jc: nonempty source cell index in level l
        keyj = sourceIndex2Key[l][jc];
        getIX(jX, keyj);                      // jc -> keyj -> jX
        if (isNeighbor(iX/2, jX/2)) {         // if ic,jc's parents are neighbor
          if (!isNeighbor(iX, jX)) {          // and ic,jc are non-neighbor
            getX(Xi, iX, l);
            getX(Xj, jX, l);
            kernel.M2L(Xi, L[keyi+keyLevelOffset[l]], Xj, M[keyj+keyLevelOffset[l]]);
          }         
        }
      }
    }
  }
  
  // L2L
  for (l=3; l<=level; l++) {                  // l: current level
    for (ic=0; ic<targetNonEmpty[l]; ic++) {  // ic: nonempty target cell index in level l
      key = targetIndex2Key[l][ic];           // key: i's Morton key
      getX(Xi, key, l);                       // Xi: child cell center
      getX(XI, key/8, l-1);                   // XI: parent cell center
      kernel.L2L(Xi, L[key+keyLevelOffset[l]], XI, L[key/8+keyLevelOffset[l-1]]);
    }
  }

  // L2P
  for (ic=0; ic<targetNonEmpty[level]; ic++) {
    key = targetIndex2Key[level][ic];
    getX(Xi, key, level);
    Ti = &targets[targetOffset[ic]];
    ni = targetOffset[ic+1] - targetOffset[ic];
    kernel.L2P(Ti, ni, Xi, L[key+keyLevelOffset[level]]);
  }

  // P2P
  for (ic=0; ic<targetNonEmpty[level]; ic++) {
    getIX(iX, targetIndex2Key[level][ic]);
    for (jc=0; jc<sourceNonEmpty[level]; jc++) {
      getIX(jX, sourceIndex2Key[level][jc]);
      if (isNeighbor(iX, jX)) {
        Ti = &targets[targetOffset[ic]];
        Sj = &sources[sourceOffset[jc]];
        ni = targetOffset[ic+1] - targetOffset[ic];
        nj = sourceOffset[jc+1] - sourceOffset[jc];
        kernel.P2P(Ti, ni, Sj, nj);
      }
    }
  }
  
  // verify
  Ti = &targets2[0];
  Sj = &sources[0];
  kernel.P2P(Ti, N, Sj, N);

  std::fstream file;
  file.open("kernel.dat", std::ios::out | std::ios::app);
  double potDif = verify.getDifScalar(targets, targets2);
  double potNrm = verify.getNrmScalar(targets);
  double accDif = verify.getDifVector(targets, targets2);
  double accNrm = verify.getNrmVector(targets);
  std::cout << args.P << " " << std::sqrt(potDif/potNrm) << "  " << std::sqrt(accDif/accNrm) << std::endl;
  double potRel = std::sqrt(potDif/potNrm);
  double accRel = std::sqrt(accDif/accNrm);
  verify.print("Rel. L2 Error (pot)",potRel);
  verify.print("Rel. L2 Error (acc)",accRel);
  file << args.P << " " << std::sqrt(potDif/potNrm) << "  " << std::sqrt(accDif/accNrm) << std::endl;
  file.close();
  return 0;
}
