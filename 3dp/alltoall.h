#ifndef alltoall_h
#define alltoall_h
#include "exafmm.h"

namespace exafmm {
  //! Use alltoall to get recv count and calculate displacement from it
  void getCountAndDispl(std::vector<int> & sendCount, std::vector<int> & sendDispl,
                        std::vector<int> & recvCount, std::vector<int> & recvDispl) {
    MPI_Alltoall(&sendCount[0], 1, MPI_INT, &recvCount[0], 1, MPI_INT, MPI_COMM_WORLD);
    for (int irank=0; irank<MPISIZE-1; irank++) {
      sendDispl[irank+1] = sendDispl[irank] + sendCount[irank];
      recvDispl[irank+1] = recvDispl[irank] + recvCount[irank];
    }
  }

  //! Alltoallv for bodies
  void alltoallBodies(Bodies & sendBodies, std::vector<int> & sendBodyCount, std::vector<int> & sendBodyDispl,
                      Bodies & recvBodies, std::vector<int> & recvBodyCount, std::vector<int> & recvBodyDispl) {
    MPI_Datatype MPI_BODY;
    MPI_Type_contiguous(sizeof(sendBodies[0]), MPI_CHAR, &MPI_BODY);
    MPI_Type_commit(&MPI_BODY);
    recvBodies.resize(recvBodyDispl[MPISIZE-1]+recvBodyCount[MPISIZE-1]);
    MPI_Alltoallv(&sendBodies[0], &sendBodyCount[0], &sendBodyDispl[0], MPI_BODY,
                  &recvBodies[0], &recvBodyCount[0], &recvBodyDispl[0], MPI_BODY, MPI_COMM_WORLD);
  }

  //! Alltoallv for cells
  void alltoallCells(Cells & sendCells, std::vector<int> & sendCellCount, std::vector<int> & sendCellDispl,
                     Cells & recvCells,  std::vector<int> & recvCellCount, std::vector<int> & recvCellDispl) {
    //! Copy cells to cell bases, cell data
    std::vector<CellBase> sendCellBases(sendCellDispl[MPISIZE-1]+sendCellCount[MPISIZE-1]);
    std::vector<complex_t> sendCellData(sendCellBases.size()*NTERM);
    for (int irank=0, ic=0, ib=0; irank<MPISIZE; irank++) {
      for (int i=sendCellDispl[irank]; i<sendCellDispl[irank]+sendCellCount[irank]; i++) {
        sendCellBases[i] = sendCells[i];
        for (int n=0; n<NTERM; n++) {
          sendCellData[ic++] = sendCells[i].M[n];
        }
      }
    }
    //! Send cell bases
    MPI_Datatype MPI_CELL_BASE;
    MPI_Type_contiguous(sizeof(sendCellBases[0]), MPI_CHAR, &MPI_CELL_BASE);
    MPI_Type_commit(&MPI_CELL_BASE);
    std::vector<CellBase> recvCellBases(recvCellDispl[MPISIZE-1]+recvCellCount[MPISIZE-1]);
    MPI_Alltoallv(&sendCellBases[0], &sendCellCount[0], &sendCellDispl[0], MPI_CELL_BASE,
                  &recvCellBases[0], &recvCellCount[0], &recvCellDispl[0], MPI_CELL_BASE, MPI_COMM_WORLD);
    //! Send cell data
    MPI_Datatype MPI_CELL_DATA;
    MPI_Type_contiguous(sizeof(complex_t)*NTERM, MPI_CHAR, &MPI_CELL_DATA);
    MPI_Type_commit(&MPI_CELL_DATA);
    std::vector<complex_t> recvCellData(recvCellBases.size()*NTERM);
    MPI_Alltoallv(&sendCellData[0], &sendCellCount[0], &sendCellDispl[0], MPI_CELL_DATA,
                  &recvCellData[0], &recvCellCount[0], &recvCellDispl[0], MPI_CELL_DATA, MPI_COMM_WORLD);
    //! Copy cell bases, cell data to cells
    recvCells.resize(recvCellBases.size());
    for (int irank=0, ic=0; irank<MPISIZE; irank++) {
      for (int i=recvCellDispl[irank]; i<recvCellDispl[irank]+recvCellCount[irank]; i++) {
        recvCells[i] = recvCellBases[i];
        recvCells[i].M.resize(NTERM, 0);
        recvCells[i].L.resize(NTERM, 0);
        for (int n=0; n<NTERM; n++) {
          recvCells[i].M[n] += recvCellData[ic++];
        }
      }
    }
  }
}
#endif
