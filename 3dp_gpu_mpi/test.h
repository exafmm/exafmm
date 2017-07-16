#ifndef test_h
#define test_h
#include "exafmm.h"

namespace exafmm {
  void P2M(Cell * C) {
    for (Body * B=C->body; B!=C->body+C->numBodies; ++B) C->M[0] += B->q;
  }

  void M2M(Cell * Ci) {
    for (Cell * Cj=Ci->child; Cj!=Ci->child+Ci->numChilds; ++Cj) Ci->M[0] += Cj->M[0];
  }

  inline void M2L(Cell * Ci, Cell * Cj) {
    Ci->L[0] += Cj->M[0];
  }

  void L2L(Cell * Cj) {
    for (Cell * Ci=Cj->child; Ci!=Cj->child+Cj->numChilds; ++Ci) Ci->L[0] += Cj->L[0];
  }

  void L2P(Cell * C) {
    for (Body * B=C->body; B!=C->body+C->numBodies; ++B) B->p += std::real(C->L[0]);
  }

  void P2P(Cell * Ci, Cell * Cj) {
    Body * Bi = Ci->body;
    Body * Bj = Cj->body;
    int ni = Ci->numBodies;
    int nj = Cj->numBodies;
    for (int i=0; i<ni; i++) {
      for (int j=0; j<nj; j++) {
        Bi[i].p += Bj[j].q;
      }
    }
  }

  void upwardPass(Cell * Ci) {
    for (Cell * Cj=Ci->child; Cj!=Ci->child+Ci->numChilds; ++Cj) {
      upwardPass(Cj);
    }
    Ci->M.resize(1, 0.0);
    Ci->L.resize(1, 0.0);
    if (Ci->numChilds==0) P2M(Ci);
    M2M(Ci);
  }

  void horizontalPass(Cell * Ci, Cell * Cj) {
    vec3 dX;
    for (int d=0; d<3; d++) dX[d] = Ci->X[d] - Cj->X[d] - IX[d] * CYCLE;
    real_t R2 = norm(dX) * THETA * THETA;
    if (R2 > (Ci->R + Cj->R) * (Ci->R + Cj->R)) {
      M2L(Ci, Cj);
    } else if (Ci->numChilds == 0 && Cj->numChilds == 0) {
      P2P(Ci, Cj);
    } else if (Cj->numChilds == 0 || (Ci->R >= Cj->R && Ci->numChilds != 0)) {
      for (Cell * ci=Ci->child; ci!=Ci->child+Ci->numChilds; ci++) {
        horizontalPass(ci, Cj);
      }
    } else {
      for (Cell * cj=Cj->child; cj!=Cj->child+Cj->numChilds; cj++) {
        horizontalPass(Ci, cj);
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

  void downwardPass(Cell * Cj) {
    L2L(Cj);
    if (Cj->numChilds==0) L2P(Cj);
    for (Cell * Ci=Cj->child; Ci!=Cj->child+Cj->numChilds; ++Ci) {
      downwardPass(Ci);
    }
  }
}
#endif
