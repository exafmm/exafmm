#ifndef local_essential_tree_h
#define local_essential_tree_h
#include "alltoall.h"
#include "hilbert.h"
#include <map>
#include "timer.h"
#define SEND_ALL 0 //! Set to 1 for debugging

namespace exafmm {
  typedef std::multimap<uint64_t, Body> BodyMap;
  typedef std::map<uint64_t, Cell> CellMap;
  int LEVEL;                                    //!< Octree level used for partitioning
  std::vector<int> OFFSET;                      //!< Offset of Hilbert index for partitions

  //! Distance between cell center and edge of a remote domain
  real_t getDistance(Cell * C, int irank) {
    real_t distance = R0;
    real_t R = R0 / (1 << LEVEL);
    for (int key=OFFSET[irank]; key<OFFSET[irank+1]; key++) {
      ivec3 iX = get3DIndex(key, LEVEL);
      vec3 X = getCoordinates(iX, LEVEL);
      vec3 Xmin = X - R;
      vec3 Xmax = X + R;
      vec3 dX;
      if (IMAGES == 0) {
        for (int d=0; d<3; d++) {
          dX[d] = (C->X[d] > Xmax[d]) * (C->X[d] - Xmax[d]) + (C->X[d] < Xmin[d]) * (C->X[d] - Xmin[d]);
        }
        distance = std::min(distance, norm(dX));
      } else {
        for (iX[0]=-1; iX[0]<=1; iX[0]++) {
          for (iX[1]=-1; iX[1]<=1; iX[1]++) {
            for (iX[2]=-1; iX[2]<=1; iX[2]++) {
              for (int d=0; d<3; d++) {
                dX[d] = (C->X[d] + iX[d] * CYCLE > Xmax[d]) * (C->X[d] + iX[d] * CYCLE - Xmax[d]) +
                  (C->X[d] + iX[d] * CYCLE < Xmin[d]) * (C->X[d] + iX[d] * CYCLE - Xmin[d]);
              }
              distance = std::min(distance, norm(dX));
            }
          }
        }
      }
    }
    return distance;
  }

  //! Recursive call to pre-order tree traversal for selecting cells to send
  void selectCells(Cell * Cj, int irank, Bodies & bodyBuffer, std::vector<int> & sendBodyCount,
                   Cells & cellBuffer, std::vector<int> & sendCellCount) {
    real_t R = getDistance(Cj, irank);
    real_t R2 = R * R * THETA * THETA;
    sendCellCount[irank]++;
    cellBuffer.push_back(*Cj);
    if (R2 <= (Cj->R + Cj->R) * (Cj->R + Cj->R)) {
      if (Cj->numChilds == 0) {
        sendBodyCount[irank] += Cj->numBodies;
        for (int b=0; b<Cj->numBodies; b++) {
          bodyBuffer.push_back(Cj->body[b]);
        }
      } else {
        for (Cell * Ci=Cj->child; Ci!=Cj->child+Cj->numChilds; Ci++) {
          selectCells(Ci, irank, bodyBuffer, sendBodyCount, cellBuffer, sendCellCount);
        }
      }
    }
  }

  void whatToSend(Cells & cells, Bodies & bodyBuffer, std::vector<int> & sendBodyCount,
                  Cells & cellBuffer, std::vector<int> & sendCellCount) {
#if SEND_ALL //! Send everything (for debugging)
    for (int irank=0; irank<MPISIZE; irank++) {
      sendCellCount[irank] = cells.size();
      for (size_t i=0; i<cells.size(); i++) {
        if (cells[i].numChilds == 0) {
          sendBodyCount[irank] += cells[i].numBodies;
          for (int b=0; b<cells[i].numBodies; b++) {
            bodyBuffer.push_back(cells[i].body[b]);
          }
        }
      }
      cellBuffer.insert(cellBuffer.end(), cells.begin(), cells.end());
    }
#else //! Send only necessary cells
    for (int irank=0; irank<MPISIZE; irank++) {
      selectCells(&cells[0], irank, bodyBuffer, sendBodyCount, cellBuffer, sendCellCount);
    }
#endif
  }

  //! Reapply Ncrit recursively to account for bodies from other ranks
  void reapplyNcrit(BodyMap & bodyMap, CellMap & cellMap, uint64_t key) {
    bool noChildSent = true;
    for (int i=0; i<8; i++) {
      uint64_t childKey = getChild(key) + i;
      if (cellMap.find(childKey) != cellMap.end()) noChildSent = false;
    }
    if (cellMap[key].numBodies <= NCRIT && noChildSent) {
      cellMap[key].numChilds = 0;
      cellMap[key].numBodies = bodyMap.count(key);
      return;
    }
    int level = getLevel(key);
    int counter[8] = {0};
    //! Assign key of child to bodyMap
    std::pair<BodyMap::iterator,BodyMap::iterator> range = bodyMap.equal_range(key);
    Bodies bodies(bodyMap.count(key));
    size_t b = 0;
    for (BodyMap::iterator B=range.first; B!=range.second; B++, b++) {
      bodies[b] = B->second;
    }
    for (b=0; b<bodies.size(); b++) {
      ivec3 iX = get3DIndex(bodies[b].X, level+1);
      uint64_t childKey = getKey(iX, level+1);
      int octant = getOctant(childKey);
      counter[octant]++;
      bodies[b].key = childKey;
      bodyMap.insert(std::pair<uint64_t, Body>(childKey, bodies[b]));
    }
    if (bodyMap.count(key) != 0) bodyMap.erase(key);
    //! Create a new cell if it didn't exist
    for (int i=0; i<8; i++) {
      uint64_t childKey = getChild(key) + i;
      if (counter[i] != 0) {
        if (cellMap.find(childKey) == cellMap.end()) {
          Cell cell;
          cell.numBodies = counter[i];
          cell.numChilds = 0;
          ivec3 iX = get3DIndex(childKey);
          cell.X = getCoordinates(iX, level+1);
          cell.R = R0 / (1 << (level+1));
          cell.key = childKey;
          cell.M.resize(NTERM, 0.0);
          cell.L.resize(NTERM, 0.0);
          cellMap[childKey] = cell;
        } else {
          cellMap[childKey].numBodies += counter[i];
          for (int n=0; n<NTERM; n++) cellMap[childKey].M[n] = 0;
        }
      }
      if (cellMap.find(childKey) != cellMap.end()) {
        reapplyNcrit(bodyMap, cellMap, childKey);
      }
    }
    //! Update number of bodies and child cells
    int numBodies = 0;
    int numChilds = 0;
    for (int i=0; i<8; i++) {
      uint64_t childKey = getChild(key) + i;
      if (cellMap.find(childKey) != cellMap.end()) {
        numBodies += cellMap[childKey].numBodies;
        numChilds++;
      }
    }
    if (numChilds == 0) numBodies = bodyMap.count(key);
    cellMap[key].numBodies = numBodies;
    cellMap[key].numChilds = numChilds;
  }

  //! Check integrity of local essential tree
  void sanityCheck(BodyMap & bodyMap, CellMap & cellMap, uint64_t key) {
    Cell cell = cellMap[key];
    assert(cell.key == key);
    if (cell.numChilds == 0) assert(cell.numBodies == int(bodyMap.count(key)));
    if (bodyMap.count(key) != 0) {
      assert(cell.numChilds == 0);
      std::pair<BodyMap::iterator,BodyMap::iterator> range = bodyMap.equal_range(key);
      for (BodyMap::iterator B=range.first; B!=range.second; B++) {
        assert(B->second.key == key);
      }
    }
    int numBodies = 0;
    int numChilds = 0;
    for (int i=0; i<8; i++) {
      uint64_t childKey = getChild(key) + i;
      if (cellMap.find(childKey) != cellMap.end()) {
        sanityCheck(bodyMap, cellMap, childKey);
        numBodies += cellMap[childKey].numBodies;
        numChilds++;
      }
    }
    assert((cell.numBodies == numBodies) || (numBodies == 0));
    assert((cell.numChilds == numChilds));
  }

  //! Build cells of LET recursively
  void buildCells(BodyMap & bodyMap, CellMap & cellMap, uint64_t key, Bodies & bodies, Cell * cell, Cells & cells) {
    *cell = cellMap[key];
    if (bodyMap.count(key) != 0) {
      std::pair<BodyMap::iterator,BodyMap::iterator> range = bodyMap.equal_range(key);
      bodies.resize(bodies.size()+cell->numBodies);
      Body * body = &bodies.back() - cell->numBodies + 1;
      cell->body = body;
      int b = 0;
      for (BodyMap::iterator B=range.first; B!=range.second; B++, b++) {
        body[b] = B->second;
      }
    } else {
      cell->body = NULL;
    }
    if (cell->numChilds != 0) {
      cells.resize(cells.size()+cell->numChilds);
      Cell * child = &cells.back() - cell->numChilds + 1;
      cell->child = child;
      for (int i=0, c=0; i<8; i++) {
        uint64_t childKey = getChild(key) + i;
        if (cellMap.find(childKey) != cellMap.end()) {
          buildCells(bodyMap, cellMap, childKey, bodies, &child[c++], cells);
        }
      }
    } else {
      cell->child = NULL;
    }
    if (cell->numChilds != 0) cell->body = cell->child->body;
  }

  //! Build local essential tree
  void buildLocalEssentialTree(Bodies & recvBodies, Cells & recvCells, Bodies & bodies, Cells & cells) {
    BodyMap bodyMap;
    CellMap cellMap;
    //! Insert bodies to multimap
    for (size_t i=0; i<recvBodies.size(); i++) {
      bodyMap.insert(std::pair<uint64_t, Body>(recvBodies[i].key, recvBodies[i]));
    }
    //! Insert cells to map and merge cells
    for (size_t i=0; i<recvCells.size(); i++) {
      uint64_t key = recvCells[i].key;
      if (cellMap.find(key) == cellMap.end()) {
        cellMap[key] = recvCells[i];
      } else {
        for (int n=0; n<NTERM; n++) {
          cellMap[key].M[n] += recvCells[i].M[n];
        }
        cellMap[key].numBodies += recvCells[i].numBodies;
      }
    }
    //! Reapply Ncrit recursively to account for bodies from other ranks
    reapplyNcrit(bodyMap, cellMap, 0);
    //! Check integrity of local essential tree
    sanityCheck(bodyMap, cellMap, 0);
    //! Copy bodyMap to bodies
    bodies.clear();
    bodies.reserve(bodyMap.size());
    //! Build cells of LET recursively
    cells.reserve(cellMap.size());
    cells.resize(1);
    buildCells(bodyMap, cellMap, 0, bodies, &cells[0], cells);
    //! Check correspondence between vector and map sizes
    assert(bodies.size() == bodyMap.size());
    assert(cells.size() == cellMap.size());
  }

  //! MPI communication for local essential tree
  void localEssentialTree(Bodies & bodies, Cells & cells) {
    std::vector<int> sendBodyCount(MPISIZE, 0);
    std::vector<int> recvBodyCount(MPISIZE, 0);
    std::vector<int> sendBodyDispl(MPISIZE, 0);
    std::vector<int> recvBodyDispl(MPISIZE, 0);
    std::vector<int> sendCellCount(MPISIZE, 0);
    std::vector<int> recvCellCount(MPISIZE, 0);
    std::vector<int> sendCellDispl(MPISIZE, 0);
    std::vector<int> recvCellDispl(MPISIZE, 0);
    Bodies sendBodies, recvBodies;
    Cells sendCells, recvCells;
    //! Decide which cells & bodies to send
    whatToSend(cells, sendBodies, sendBodyCount, sendCells, sendCellCount);
    //! Use alltoall to get recv count and calculate displacement (defined in alltoall.h)
    getCountAndDispl(sendBodyCount, sendBodyDispl, recvBodyCount, recvBodyDispl);
    getCountAndDispl(sendCellCount, sendCellDispl, recvCellCount, recvCellDispl);
    //! Alltoallv for cells (defined in alltoall.h)
    alltoallCells(sendCells, sendCellCount, sendCellDispl, recvCells, recvCellCount, recvCellDispl);
    //! Alltoallv for bodies (defined in alltoall.h)
    alltoallBodies(sendBodies, sendBodyCount, sendBodyDispl, recvBodies, recvBodyCount, recvBodyDispl);
    //! Build local essential tree
    buildLocalEssentialTree(recvBodies, recvCells, bodies, cells);
  }
}
#endif
