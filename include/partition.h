#ifndef partition_h
#define partition_h
#include "alltoall.h"
#include "hilbert.h"

namespace exafmm {
  //! Allreduce local bounds to get global bounds
  void allreduceBounds(Bodies & bodies) {
    vec3 localXmin, localXmax, globalXmin, globalXmax;
    localXmin = bodies[0].X;
    localXmax = bodies[0].X;
    for (size_t b=0; b<bodies.size(); b++) {
      localXmin = min(bodies[b].X, localXmin);
      localXmax = max(bodies[b].X, localXmax);
    }
    MPI_Allreduce(&localXmin[0], &globalXmin[0], 3, MPI_REAL_T, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&localXmax[0], &globalXmax[0], 3, MPI_REAL_T, MPI_MAX, MPI_COMM_WORLD);
    X0 = (globalXmax + globalXmin) / 2;
    R0 = max(X0-globalXmin);
    R0 *= 1.00001;
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

  void partition(Bodies & bodies) {
    allreduceBounds(bodies);
    const int numBodies = bodies.size();
    const int numBins = 1 << 3 * LEVEL;
    std::vector<int> localHist(numBins, 0);
    std::vector<int> key(numBodies);
    std::vector<int> index(numBodies);
    //! Get local histogram of hilbert key bins
    for (int b=0; b<numBodies; b++) {
      ivec3 iX = get3DIndex(bodies[b].X, LEVEL);
      key[b] = getKey(iX, LEVEL, false);
      index[b] = b;
      localHist[key[b]]++;
    }
    //! Sort bodies according to keys
    std::vector<int> key2 = key;
    radixsort(key, index, numBodies);
    Bodies buffer = bodies;
    for (int b=0; b<numBodies; b++) {
      bodies[b] = buffer[index[b]];
      key[b] = key2[index[b]];
    }
    //! Get Global histogram of hilbert key bins
    std::vector<int> globalHist(numBins);
    MPI_Allreduce(&localHist[0], &globalHist[0], numBins, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    //! Calculate offset of global histogram for each rank
    OFFSET.resize(MPISIZE+1);
    OFFSET[0] = 0;
    for (int i=0, irank=0, count=0; i<numBins; i++) {
      count += globalHist[i];
      if (irank * numBodies < count) {
        OFFSET[irank] = i;
        irank++;
      }
    }
    OFFSET[MPISIZE] = numBins;
    std::vector<int> sendBodyCount(MPISIZE, 0);
    std::vector<int> recvBodyCount(MPISIZE, 0);
    std::vector<int> sendBodyDispl(MPISIZE, 0);
    std::vector<int> recvBodyDispl(MPISIZE, 0);
    //! Use the offset as the splitter for partitioning
    for (int irank=0, b=0; irank<MPISIZE; irank++) {
      while (b < int(bodies.size()) && key[b] < OFFSET[irank+1]) {
        sendBodyCount[irank]++;
        b++;
      }
    }
    //! Use alltoall to get recv count and calculate displacement from it
    getCountAndDispl(sendBodyCount, sendBodyDispl, recvBodyCount, recvBodyDispl);
    buffer.resize(recvBodyDispl[MPISIZE-1]+recvBodyCount[MPISIZE-1]);
    //! Alltoallv for bodies (defined in alltoall.h)
    alltoallBodies(bodies, sendBodyCount, sendBodyDispl, buffer, recvBodyCount, recvBodyDispl);
    bodies = buffer;
  }

  //! Shift bodies among MPI rank round robin
  void shiftBodies(Bodies & bodies) {
    int newSize;
    int oldSize = bodies.size();
    const int isend = (MPIRANK + 1          ) % MPISIZE;
    const int irecv = (MPIRANK - 1 + MPISIZE) % MPISIZE;
    MPI_Request sreq,rreq;
    MPI_Datatype MPI_BODY;
    MPI_Type_contiguous(sizeof(bodies[0]), MPI_CHAR, &MPI_BODY);
    MPI_Type_commit(&MPI_BODY);
    MPI_Isend(&oldSize, 1, MPI_INT, irecv, 0, MPI_COMM_WORLD, &sreq);
    MPI_Irecv(&newSize, 1, MPI_INT, isend, 0, MPI_COMM_WORLD, &rreq);
    MPI_Wait(&sreq, MPI_STATUS_IGNORE);
    MPI_Wait(&rreq, MPI_STATUS_IGNORE);
    Bodies buffer = bodies;
    bodies.resize(newSize);
    MPI_Isend(&buffer[0], oldSize, MPI_BODY, irecv, 1, MPI_COMM_WORLD, &sreq);
    MPI_Irecv(&bodies[0], newSize, MPI_BODY, isend, 1, MPI_COMM_WORLD, &rreq);
    MPI_Wait(&sreq, MPI_STATUS_IGNORE);
    MPI_Wait(&rreq, MPI_STATUS_IGNORE);
  }

  //! Allgather bodies
  void gatherBodies(Bodies & bodies) {
    std::vector<int> recvCount(MPISIZE, 0);
    std::vector<int> recvDispl(MPISIZE, 0);
    int sendCount = bodies.size();
    MPI_Datatype MPI_BODY;
    MPI_Type_contiguous(sizeof(bodies[0]), MPI_CHAR, &MPI_BODY);
    MPI_Type_commit(&MPI_BODY);
    MPI_Allgather(&sendCount, 1, MPI_INT, &recvCount[0], 1, MPI_INT, MPI_COMM_WORLD);
    for (int irank=0; irank<MPISIZE-1; irank++) {
      recvDispl[irank+1] = recvDispl[irank] + recvCount[irank];
    }
    Bodies buffer = bodies;
    bodies.resize(recvDispl[MPISIZE-1]+recvCount[MPISIZE-1]);
    MPI_Allgatherv(&buffer[0], sendCount, MPI_BODY,
                   &bodies[0], &recvCount[0], &recvDispl[0], MPI_BODY, MPI_COMM_WORLD);
  }
}
#endif
