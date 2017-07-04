#ifndef partition_h
#define partition_h
#include "base_mpi.h"
#include "exafmm.h"

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
      numLevels = 0;
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
      buffer = bodies;                                          // Resize sort buffer
      for (int level=0; level<numLevels; level++) {             // Loop over levels of MPI rank binary tree
        int numPartitions = rankColor[mpisize-1] + 1;           //  Number of partitions in current level
        for (int ipart=0; ipart<numPartitions; ipart++) {       //  Loop over partitions
          int irank = rankMap[ipart];                           //   Map partition to MPI rank
          Bounds bounds = rankBounds[irank];                    //   Bounds of current rank
          int direction = 0;                                    //   Initialize direction of partitioning
          real_t length = 0;                                    //   Initialize length of partition
          for (int d=0; d<3; d++) {                             //   Loop over dimensions
            if (length < (bounds.Xmax[d] - bounds.Xmin[d])) {   //    If present dimension is longer
              direction = d;                                    //     Update direction to current dimension
              length = (bounds.Xmax[d] - bounds.Xmin[d]);       //     Update length to current dimension
            }                                                   //    End if for longer dimension
          }                                                     //   End loop over dimensions
          int rankSplit = rankCount[irank] / 2;                 //   MPI rank splitter
          int oldRankCount = rankCount[irank];                  //   Old MPI rank count
          int bodyBegin = 0;                                    //   Initialize body begin index
          int bodyEnd = sendCount[irank];                       //   Initialize body end index
          B_iter B = bodies.begin() + sendDispl[irank];         //   Body begin iterator of current partition
          float localWeightSum = 0;                             //   Initialize local sum of weights in current partition
          for (int b=bodyBegin; b<bodyEnd; b++) {               //    Loop over bodies in current partition
            localWeightSum += B[b].weight;                      //     Add weights of body to local sum
          }                                                     //    End loop over bodies in current partition
          float globalWeightSum;                                //   Declare global sum of weights in current partition
          MPI_Allreduce(&localWeightSum, &globalWeightSum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);// Reduce sum of weights
          float globalSplit = globalWeightSum * rankSplit / oldRankCount;// Global weight splitter index
          float globalOffset = 0;                               //   Initialize global weight offset
          real_t xmax = bounds.Xmax[direction];                 //   Upper bound of partition
          real_t xmin = bounds.Xmin[direction];                 //   Lower bound of partition
          real_t dx = (xmax - xmin) / numBins;                  //   Length of bins
          if (rankSplit > 0) {                                  //   If the partition requires splitting
            for (int binRefine=0; binRefine<3; binRefine++) {   //    Loop for bin refinement
              for (int ibin=0; ibin<numBins; ibin++) {          //     Loop over bins
                countHist[ibin] = 0;                            //      Initialize body count histogram
                weightHist[ibin] = 0;                           //      Initialize body weight histogram
              }                                                 //     End loop over bins
              for (int b=bodyBegin; b<bodyEnd; b++) {           //     Loop over bodies
                real_t x = B[b].X[direction];                   //      Coordinate of body in current direction
                int ibin = (x - xmin + EPS) / (dx + EPS);       //      Assign bin index to body
                countHist[ibin]++;                              //      Increment body count histogram
                weightHist[ibin] += B[b].weight;                //      Increment body weight histogram
              }                                                 //     End loop over bodies
              scanHist[0] = countHist[0];                       //     Initialize scan array
              for (int ibin=1; ibin<numBins; ibin++) {          //     Loop over bins
                scanHist[ibin] = scanHist[ibin-1] + countHist[ibin]; //   Inclusive scan for body count histogram
              }                                                 //     End loop over bins
              for (int b=bodyEnd-1; b>=bodyBegin; b--) {        //     Loop over bodies backwards
                real_t x = B[b].X[direction];                   //      Coordinate of body in current direction
                int ibin = (x - xmin + EPS) / (dx + EPS);       //      Assign bin index to body
                scanHist[ibin]--;                               //      Decrement scanned histogram
                int bnew = scanHist[ibin] + bodyBegin;          //      New index of sorted body
                buffer[bnew] = B[b];                            //      Copy body to sort buffer
              }                                                 //     End loop over bodies backwards
              for (int b=bodyBegin; b<bodyEnd; b++) {           //     Loop over bodies
                B[b] = buffer[b];                               //      Copy back bodies from buffer
              }                                                 //     End loop over bodies
              MPI_Allreduce(&weightHist[0], &globalHist[0], numBins, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);// Reduce weight histogram
              int splitBin = 0;                                 //     Initialize bin splitter
              while (globalOffset < globalSplit) {              //     While scan of global histogram is less than splitter
                globalOffset += globalHist[splitBin];           //      Scan global histogram
                splitBin++;                                     //      Increment bin count
              }                                                 //     End while for global histogram scan
              splitBin--;                                       //     Move back one bin
              globalOffset -= globalHist[splitBin];             //     Decrement offset accordingly
              xmax = xmin + (splitBin + 1) * dx;                //     Zoom in to current bin by redefining upper
              xmin = xmin + splitBin * dx;                      //     and lower bounds of partition
              dx = (xmax - xmin) / numBins;                     //     Update length of partition accordingly
              bodyBegin += scanHist[splitBin];                  //     Update body begin index
              bodyEnd = bodyBegin + countHist[splitBin];        //     Update body end index
            }                                                   //    End loop for bin refinement
          }                                                     //   End if for splitting partition
          int rankBegin = rankDispl[irank];                     //   Save current range of MPI ranks
          int rankEnd = rankBegin + rankCount[irank];           //   so that they don't get overwritten
          for (irank=rankBegin; irank<rankEnd; irank++) {       //   Loop over current range of MPI ranks
            rankSplit = rankCount[irank] / 2;                   //    MPI rank splitter
            if (irank - rankDispl[irank] < rankSplit) {         //    If on left side of splitter
              rankCount[irank] = rankSplit;                     //     Count is the splitter index
              rankColor[irank] = rankColor[irank] * 2;          //     Color is doubled
              rankBounds[irank].Xmax[direction] = xmin;         //     Update right bound with splitter
              sendCount[irank] = bodyBegin;                     //     Send body count is the begin index
            } else {                                            //    If on right side of splitter
              rankDispl[irank] += rankSplit;                    //     Update displacement with splitter index
              rankCount[irank] -= rankSplit;                    //     Count is remainder of splitter index
              rankColor[irank] = rankColor[irank] * 2 + 1;      //     Color is doubled plus one
              rankBounds[irank].Xmin[direction] = xmin;         //     Update left bound
              sendDispl[irank] += bodyBegin;                    //     Increment send body displacement
              sendCount[irank] -= bodyBegin;                    //     Decrement send body count
            }                                                   //    End if for side of splitter
            if (level == numLevels-1) rankColor[irank] = rankDispl[irank]; // Special case for final rank color
            rankKey[irank] = irank - rankDispl[irank];          //    Rank key is determined from rank displacement
          }                                                     //   End loop over current range of MPI ranks
        }                                                       //  End loop over partitions
        int ipart = 0;                                          //  Initialize partition index
        for (int irank=0; irank<mpisize; irank++) {             //  Loop over MPI ranks
          if (rankKey[irank] == 0) {                            //   If rank key is zero it's a new group
            rankMap[ipart] = rankDispl[irank];                  //    Update rank map with rank displacement
            ipart++;                                            //    Increment partition index
          }                                                     //   End if for new group
        }                                                       //  End loop over MPI ranks
      }                                                         // End loop over levels of MPI rank binary tree
      B_iter B = bodies.begin();                                // Body begin iterator
      for (int irank=0; irank<mpisize; irank++) {               // Loop over MPI ranks
        int bodyBegin = sendDispl[irank];                       //  Body begin index for current rank
        int bodyEnd = bodyBegin + sendCount[irank];             //  Body end index for current rank
        for (int b=bodyBegin; b<bodyEnd; b++, B++) {            //  Loop over bodies in current rank
          B->irank = irank;                                     //   Copy MPI rank to body
        }                                                       //  End loop over bodies in current rank
      }                                                         // End loop over MPI ranks
      return rankBounds[mpirank];                               // Return local bounds
    }
  };
}
#endif
