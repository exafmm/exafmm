#ifndef traverse_lazy_h
#define traverse_lazy_h
#include "exafmm.h"
#include "kernel.h"

namespace exafmm {
  //! Recursive call to post-order tree traversal for upward pass
  void upwardPass(Cell * Ci) {
    for (Cell * Cj=Ci->child; Cj!=Ci->child+Ci->numChilds; Cj++) {
#pragma omp task untied if(Cj->numBodies > 100)
      upwardPass(Cj);
    }
#pragma omp taskwait
    Ci->M.resize(NTERM, 0.0);
    Ci->L.resize(NTERM, 0.0);
    if(Ci->numChilds==0) P2M(Ci);
    M2M(Ci);
  }

  //! Upward pass interface
  void upwardPass(Cells & cells) {
#pragma omp parallel
#pragma omp single nowait
    upwardPass(&cells[0]);
  }

  //! 3-D to 1-D periodic index
  int periodic1D(int * IX) {
    return IX[0] + 1 + 3 * (IX[1] + 1) + 9 * (IX[2] + 1);
  }

  //! 1-D to 3-D periodic index
  void periodic3D(int i, int * IX) {
    IX[0] = (i % 3) - 1;
    IX[1] = ((i / 3) % 3) - 1;
    IX[2] = (i / 9) - 1;
  }

  //! Recursive call to dual tree traversal for list construction
  void getList(Cell * Ci, Cell * Cj) {
    vec3 dX;
    for (int d=0; d<3; d++) dX[d] = Ci->X[d] - Cj->X[d] - IX[d] * CYCLE;
    real_t R2 = norm(dX) * THETA * THETA;
    if (R2 > (Ci->R + Cj->R) * (Ci->R + Cj->R)) {
      Ci->listM2L.push_back(Cj);
      Ci->periodicM2L.push_back(periodic1D(IX));
    } else if (Ci->numChilds == 0 && Cj->numChilds == 0) {
      Ci->listP2P.push_back(Cj);
      Ci->periodicP2P.push_back(periodic1D(IX));
    } else if (Cj->numChilds == 0 || (Ci->R >= Cj->R && Ci->numChilds != 0)) {
      for (Cell * ci=Ci->child; ci!=Ci->child+Ci->numChilds; ci++) {
        getList(ci, Cj);
      }
    } else {
      for (Cell * cj=Cj->child; cj!=Cj->child+Cj->numChilds; cj++) {
        getList(Ci, cj);
      }
    }
  }

  //! Evaluate M2L, P2P kernels
  void evaluate(Cells & cells) {
#pragma omp parallel for
    for (size_t i=0; i<cells.size(); i++) {
      for (size_t j=0; j<cells[i].listM2L.size(); j++) {
        periodic3D(cells[i].periodicM2L[j],IX);
	M2L(&cells[i],cells[i].listM2L[j]);
      }
      for (size_t j=0; j<cells[i].listP2P.size(); j++) {
        periodic3D(cells[i].periodicP2P[j],IX);
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
    Ci->child = &pcells[0];
    Ci->numChilds = 26;
    for (int level=0; level<IMAGES-1; level++) {
      for (int ix=-1; ix<=1; ix++) {
        for (int iy=-1; iy<=1; iy++) {
          for (int iz=-1; iz<=1; iz++) {
            if (ix != 0 || iy != 0 || iz != 0) {
              for (int cx=-1; cx<=1; cx++) {
                for (int cy=-1; cy<=1; cy++) {
                  for (int cz=-1; cz<=1; cz++) {
                    IX[0] = ix * 3 + cx;
                    IX[1] = iy * 3 + cy;
                    IX[2] = iz * 3 + cz;
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
              Cj->X[0] = Ci->X[0] + ix * CYCLE;
              Cj->X[1] = Ci->X[1] + iy * CYCLE;
              Cj->X[2] = Ci->X[2] + iz * CYCLE;
              Cj->M = Ci->M;
              Cj++;
            }
          }
        }
      }
      M2M(Ci);
      CYCLE *= 3;
    }
  }

  //! Horizontal pass interface
  void horizontalPass(Cells & icells, Cells & jcells) {
    if (IMAGES == 0) {
      getList(&icells[0], &jcells[0]);
      evaluate(icells);
    } else {
      for (IX[0]=-1; IX[0]<=1; IX[0]++) {
        for (IX[1]=-1; IX[1]<=1; IX[1]++) {
          for (IX[2]=-1; IX[2]<=1; IX[2]++) {
            getList(&icells[0], &jcells[0]);
          }
        }
      }
      evaluate(icells);
      real_t saveCycle = CYCLE;
      periodic(&icells[0], &jcells[0]);
      CYCLE = saveCycle;
    }
  }

  //! Recursive call to pre-order tree traversal for downward pass
  void downwardPass(Cell * Cj) {
    L2L(Cj);
    if (Cj->numChilds==0) L2P(Cj);
    for (Cell * Ci=Cj->child; Ci!=Cj->child+Cj->numChilds; Ci++) {
#pragma omp task untied if(Ci->numBodies > 100)
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
    Ci->body = &bodies[0];
    Ci->numBodies = bodies.size();
    Cj->body = &jbodies[0];
    Cj->numBodies = jbodies.size();
    int prange = 0;
    for (int i=0; i<IMAGES; i++) {
      prange += int(powf(3.,i));
    }
#pragma omp parallel for collapse(3)
    for (int ix=-prange; ix<=prange; ix++) {
      for (int iy=-prange; iy<=prange; iy++) {
        for (int iz=-prange; iz<=prange; iz++) {
          IX[0] = ix;
          IX[1] = iy;
          IX[2] = iz;
          P2P(Ci, Cj);
        }
      }
    }
  }
}
#endif
