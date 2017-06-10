#ifndef partition_h
#define partition_h
#include "base_mpi.h"
#include "exafmm.h"
#include <vector>

namespace exafmm {
  //! Handles all the partitioning of domains
  class Partition : public BaseMPI {
  private:
    const int numBins;                                    //!< Number of sampling bins
    int numLevels;                                        //!< Levels of MPI rank binary tree
    std::vector<int> rankDispl;                           //!< Displacement of MPI rank group
    std::vector<int> rankCount;                           //!< Size of MPI rank group
    std::vector<int> rankColor;                           //!< Color of MPI rank group
    std::vector<int> rankKey;                             //!< Key of MPI rank group
    std::vector<int> rankMap;                             //!< Map partition to MPI rank group
    std::vector<int> sendDispl;                           //!< Displacement of bodies to send per rank
    std::vector<int> sendCount;                           //!< Count of bodies to send per rank
    std::vector<int> scanHist;                            //!< Scan of histogram
    std::vector<int> countHist;                           //!< Body count histogram
    std::vector<float> weightHist;                        //!< Body weight histogram
    std::vector<float> globalHist;                        //!< Global body weight histogram
    std::vector<Bounds> rankBounds;                       //!< Bounds of each rank
    Bodies buffer;                                        //!< MPI communication buffer for bodies

  public:
    //! Constructor
    Partition() : numBins(16) {
      int size = mpisize - 1;
      int numLevels = 0;
      while (size > 0) {                                  // Calculate number of levels of MPI rank binary tree
        size >>= 1;
        numLevels++;
      }
      rankDispl.resize(mpisize);
      rankCount.resize(mpisize);
      rankColor.resize(mpisize);
      rankKey.resize(mpisize);
      rankMap.resize(mpisize);
      sendDispl.resize(mpisize);
      sendCount.resize(mpisize);
      rankBounds.resize(mpisize);
      scanHist.resize(numBins);
      countHist.resize(numBins);
      weightHist.resize(numBins);
      globalHist.resize(numBins);
    } 

    //! Partitioning by orthogonal recursive bisection
    Bounds bisection(Bodies & bodies, Bounds globalBounds) {
      std::fill(rankDispl.begin(), rankDispl.end(), 0);
      std::fill(rankCount.begin(), rankCount.end(), mpisize);
      std::fill(rankColor.begin(), rankColor.end(), 0);
      std::fill(rankKey.begin(), rankKey.end(), 0);
      std::fill(rankMap.begin(), rankMap.end(), 0);
      std::fill(sendDispl.begin(), sendDispl.end(), 0);
      std::fill(sendCount.begin(), sendCount.end(), bodies.size());
      std::fill(rankBounds.begin(), rankBounds.end(), globalBounds);
    }
  };
}
#endif
