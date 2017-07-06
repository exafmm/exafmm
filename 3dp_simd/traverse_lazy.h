#ifndef traverse_lazy_h
#define traverse_lazy_h
#include "exafmm.h"
#include "kernel.h"

namespace exafmm {
  //! Recursive call to post-order tree traversal for upward pass
  void upwardPass(Cell * Ci) {
    for (Cell * Cj=Ci->CHILD; Cj!=Ci->CHILD+Ci->NCHILD; Cj++) {
#pragma omp task untied if(Cj->NBODY > 100)
      upwardPass(Cj);
    }
#pragma omp taskwait
    Ci->M.resize(NTERM, 0.0);
    Ci->L.resize(NTERM, 0.0);
    if(Ci->NCHILD==0) P2M(Ci);
    M2M(Ci);
  }

  //! Upward pass interface
  void upwardPass(Cells & cells) {
#pragma omp parallel
#pragma omp single nowait
    upwardPass(&cells[0]);
  }

  //! 3-D to 1-D periodic index
  int periodic1D(int * iX) {
    return iX[0] + 1 + 3 * (iX[1] + 1) + 9 * (iX[2] + 1);
  }

  //! 1-D to 3-D periodic index
  void periodic3D(int i, int * iX) {
    iX[0] = (i % 3) - 1;
    iX[1] = ((i / 3) % 3) - 1;
    iX[2] = (i / 9) - 1;
  }

  //! Recursive call to dual tree traversal for list construction
  void getList(Cell * Ci, Cell * Cj) {
    vec3 dX;
    for (int d=0; d<3; d++) dX[d] = Ci->X[d] - Cj->X[d] - iX[d] * cycle;
    real_t R2 = norm(dX) * theta * theta;
    if (R2 > (Ci->R + Cj->R) * (Ci->R + Cj->R)) {
      Ci->listM2L.push_back(Cj);
      Ci->periodicM2L.push_back(periodic1D(iX));
    } else if (Ci->NCHILD == 0 && Cj->NCHILD == 0) {
      Ci->listP2P.push_back(Cj);
      Ci->periodicP2P.push_back(periodic1D(iX));
    } else if (Cj->NCHILD == 0 || (Ci->R >= Cj->R && Ci->NCHILD != 0)) {
      for (Cell * ci=Ci->CHILD; ci!=Ci->CHILD+Ci->NCHILD; ci++) {
        getList(ci, Cj);
      }
    } else {
      for (Cell * cj=Cj->CHILD; cj!=Cj->CHILD+Cj->NCHILD; cj++) {
        getList(Ci, cj);
      }
    }
  }

  //! Evaluate M2L, P2P kernels
  void evaluate(Cells & cells) {
#pragma omp parallel for
    for (size_t i=0; i<cells.size(); i++) {
      for (size_t j=0; j<cells[i].listM2L.size(); j++) {
        periodic3D(cells[i].periodicM2L[j],iX);
	M2L(&cells[i],cells[i].listM2L[j]);
      }
      for (size_t j=0; j<cells[i].listP2P.size(); j++) {
        periodic3D(cells[i].periodicP2P[j],iX);
        P2P(&cells[i],cells[i].listP2P[j]);
      }
    }
  }

  //! Horizontal pass for periodic images
  void periodic(Cell * Ci0, Cell * Cj0) {
    Cells pcells(27);
    for (size_t c=0; c<pcells.size(); c++) {
      pcells[c].M.resize(NTERM, 0.0);
      pcells[c].L.resize(NTERM, 0.0);
    }
    Cell * Ci = &pcells.back();
    *Ci = *Cj0;
    Ci->CHILD = &pcells[0];
    Ci->NCHILD = 26;
    for (int level=0; level<images-1; level++) {
      for (int ix=-1; ix<=1; ix++) {
        for (int iy=-1; iy<=1; iy++) {
          for (int iz=-1; iz<=1; iz++) {
            if (ix != 0 || iy != 0 || iz != 0) {
              for (int cx=-1; cx<=1; cx++) {
                for (int cy=-1; cy<=1; cy++) {
                  for (int cz=-1; cz<=1; cz++) {
                    iX[0] = ix * 3 + cx;
                    iX[1] = iy * 3 + cy;
                    iX[2] = iz * 3 + cz;
                    M2L(Ci0, Ci);
                  }
                }
              }
            }
          }
        }
      }
      Cell * Cj = &pcells[0];
      for (int ix=-1; ix<=1; ix++) {
        for (int iy=-1; iy<=1; iy++) {
          for (int iz=-1; iz<=1; iz++) {
            if (ix != 0 || iy != 0 || iz != 0) {
              Cj->X[0] = Ci->X[0] + ix * cycle;
              Cj->X[1] = Ci->X[1] + iy * cycle;
              Cj->X[2] = Ci->X[2] + iz * cycle;
              Cj->M = Ci->M;
              Cj++;
            }
          }
        }
      }
      M2M(Ci);
      cycle *= 3;
    }
  }

  //! Horizontal pass interface
  void horizontalPass(Cells & icells, Cells & jcells) {
    if (images == 0) {
      getList(&icells[0], &jcells[0]);
      evaluate(icells);
    } else {
      for (iX[0]=-1; iX[0]<=1; iX[0]++) {
        for (iX[1]=-1; iX[1]<=1; iX[1]++) {
          for (iX[2]=-1; iX[2]<=1; iX[2]++) {
            getList(&icells[0], &jcells[0]);
          }
        }
      }
      evaluate(icells);
      real_t saveCycle = cycle;
      periodic(&icells[0], &jcells[0]);
      cycle = saveCycle;
    }
  }

  //! Recursive call to pre-order tree traversal for downward pass
  void downwardPass(Cell * Cj) {
    L2L(Cj);
    if (Cj->NCHILD==0) L2P(Cj);
    for (Cell * Ci=Cj->CHILD; Ci!=Cj->CHILD+Cj->NCHILD; Ci++) {
#pragma omp task untied if(Ci->NBODY > 100)
      downwardPass(Ci);
    }
#pragma omp taskwait
  }

  //! Downward pass interface
  void downwardPass(Cells & cells) {
#pragma omp parallel
#pragma omp single nowait
    downwardPass(&cells[0]);
  }

  //! Direct summation
  void direct(Bodies & bodies, Bodies & jbodies) {
    Cells cells(2);
    Cell * Ci = &cells[0];
    Cell * Cj = &cells[1];
    Ci->BODY = &bodies[0];
    Ci->NBODY = bodies.size();
    Cj->BODY = &jbodies[0];
    Cj->NBODY = jbodies.size();
    int prange = 0;
    for (int i=0; i<images; i++) {
      prange += int(powf(3.,i));
    }
#pragma omp parallel for collapse(3)
    for (int ix=-prange; ix<=prange; ix++) {
      for (int iy=-prange; iy<=prange; iy++) {
        for (int iz=-prange; iz<=prange; iz++) {
          iX[0] = ix;
          iX[1] = iy;
          iX[2] = iz;
          P2P(Ci, Cj);
        }
      }
    }
  }
}
#endif
