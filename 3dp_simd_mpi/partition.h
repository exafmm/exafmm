#ifndef partition_h
#define partition_h
#include "exafmm.h"
namespace exafmm {
  int LEVEL;                                    //!< Octree level used for partitioning
  std::vector<int> OFFSET;                      //!< Offset of Hilbert index for partitions

  //! Get bounding box of bodies
  void getBounds(Bodies & bodies) {
    vec3 Xmin = bodies[0].X;
    vec3 Xmax = bodies[0].X;
    for (size_t b=0; b<bodies.size(); b++) {
      Xmin = min(bodies[b].X, Xmin);
      Xmax = max(bodies[b].X, Xmax);
    }
    X0 = (Xmax + Xmin) / 2;
    R0 = fmax(max(X0-Xmin), max(Xmax-X0));
    R0 *= 1.00001;
  }

  //! Allreduce local bounds to get global bounds
  void allreduceBounds() {
    float localXmin[3], localXmax[3], globalXmin[3], globalXmax[3];
    for (int d=0; d<3; d++) {
      localXmin[d] = X0[d] - R0;
      localXmax[d] = X0[d] + R0;
    }
    MPI_Allreduce(localXmin, globalXmin, 3, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(localXmax, globalXmax, 3, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    for (int d=0; d<3; d++) {
      X0[d] = (globalXmax[d] + globalXmin[d]) / 2;
      R0 = fmax(R0, fmax(X0[d]-globalXmin[d], globalXmax[d]-X0[d]));
    }
  }

  //! Get Hilbert index from coordinate
  int getHilbert(vec3 & X, int level) {
    vec3 Xmin = X0 - R0;
    vec3 Xmax = X0 + R0;
    real_t dx = 2 * R0 / (1 << level);
    int ix = (X[0] - Xmin[0]) / dx;
    int iy = (X[1] - Xmin[1]) / dx;
    int iz = (X[2] - Xmin[2]) / dx;
    const int octantMap[8] = {0, 1, 7, 6, 3, 2, 4, 5};
    int mask = 1 << (level - 1);
    int key = 0;
    for (int i=0; i<level; i++) {
      const int bx = (ix & mask) ? 1 : 0;
      const int by = (iy & mask) ? 1 : 0;
      const int bz = (iz & mask) ? 1 : 0;
      const int octant = (bx << 2) + (by << 1) + bz;
      if(octant == 0) {
        std::swap(iy, iz);
      } else if(octant == 1 || octant == 5) {
        std::swap(ix, iy);
      } else if(octant == 4 || octant == 6) {
        ix = (ix) ^ (-1);
        iz = (iz) ^ (-1);
      } else if(octant == 3 || octant == 7) {
        ix = (ix) ^ (-1);
        iy = (iy) ^ (-1);
        std::swap(ix, iy);
      } else {
        iy = (iy) ^ (-1);
        iz = (iz) ^ (-1);
        std::swap(iy, iz);
      }
      key = (key << 3) + octantMap[octant];
      mask >>= 1;
    }
    return key;
  }

  //! Radix sort with permutation index
  void radixsort(std::vector<int> & key, std::vector<int> & value, int size) {
    const int bitStride = 8;
    const int stride = 1 << bitStride;
    const int mask = stride - 1;
    int maxKey = 0;
    int bucket[stride];
    std::vector<int> buffer(size);
    std::vector<int> permutation(size);
    for (int i=0; i<size; i++)
      if (key[i] > maxKey)
        maxKey = key[i];
    while (maxKey > 0) {
      for (int i=0; i<stride; i++)
        bucket[i] = 0;
      for (int i=0; i<size; i++)
        bucket[key[i] & mask]++;
      for (int i=1; i<stride; i++)
        bucket[i] += bucket[i-1];
      for (int i=size-1; i>=0; i--)
        permutation[i] = --bucket[key[i] & mask];
      for (int i=0; i<size; i++)
        buffer[permutation[i]] = value[i];
      for (int i=0; i<size; i++)
        value[i] = buffer[i];
      for (int i=0; i<size; i++)
        buffer[permutation[i]] = key[i];
      for (int i=0; i<size; i++)
        key[i] = buffer[i] >> bitStride;
      maxKey >>= bitStride;
    }
  }

  //! Alltoallv for partitioning bodies
  void alltoallBodies(Bodies & bodies, Bodies & buffer, std::vector<int> & key) {
    std::vector<int> sendBodyCount(MPISIZE, 0);
    for (int irank=0, b=0; irank<MPISIZE; irank++) {
      while (key[b] < OFFSET[irank+1]) {
        sendBodyCount[irank]++;
        b++;
      }
    }
    std::vector<int> recvBodyCount(MPISIZE);
    MPI_Alltoall(&sendBodyCount[0], 1, MPI_INT, &recvBodyCount[0], 1, MPI_INT, MPI_COMM_WORLD);
    std::vector<int> sendBodyDispl(MPISIZE, 0);
    std::vector<int> recvBodyDispl(MPISIZE, 0);
    for (int irank=0; irank<MPISIZE-1; irank++) {
      sendBodyDispl[irank+1] = sendBodyDispl[irank] + sendBodyCount[irank];
      recvBodyDispl[irank+1] = recvBodyDispl[irank] + recvBodyCount[irank];
    }
    MPI_Datatype MPI_BODY;
    MPI_Type_contiguous(sizeof(bodies[0])/4, MPI_INT, &MPI_BODY);
    MPI_Type_commit(&MPI_BODY);
    buffer.resize(recvBodyDispl[MPISIZE-1]+recvBodyCount[MPISIZE-1]);
    MPI_Alltoallv(&bodies[0], &sendBodyCount[0], &sendBodyDispl[0], MPI_BODY,
                  &buffer[0], &recvBodyCount[0], &recvBodyDispl[0], MPI_BODY, MPI_COMM_WORLD);
    bodies = buffer;
  }

  void partition(Bodies & bodies) {
    getBounds(bodies);
    allreduceBounds();
    const int numBodies = bodies.size();
    const int numBins = 1 << 3 * LEVEL;
    std::vector<int> localHist(numBins, 0);
    std::vector<int> key(numBodies);
    std::vector<int> index(numBodies);
    for (int b=0; b<numBodies; b++) {
      key[b] = getHilbert(bodies[b].X, LEVEL);
      index[b] = b;
      localHist[key[b]]++;
    }
    std::vector<int> key2 = key;
    radixsort(key, index, numBodies);
    Bodies buffer = bodies;
    for (int b=0; b<numBodies; b++) {
      bodies[b] = buffer[index[b]];
      key[b] = key2[index[b]];
    }
    std::vector<int> globalHist(numBins);
    MPI_Allreduce(&localHist[0], &globalHist[0], numBins, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    OFFSET.resize(MPISIZE+1);
    OFFSET[0] = 0;
    for (uint64_t i=0, irank=0, count=0; i<numBins; i++) {
      count += globalHist[i];
      if (irank * numBodies < count) {
        OFFSET[irank] = i;
        irank++;
      }
    }
    OFFSET[MPISIZE] = numBins;
    alltoallBodies(bodies, buffer, key);
  }
}
#endif
