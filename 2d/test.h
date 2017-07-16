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
    vec2 dX = Ci->X - Cj->X;
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

  void downwardPass(Cell * Cj) {
    L2L(Cj);
    if (Cj->numChilds==0) L2P(Cj);
    for (Cell * Ci=Cj->child; Ci!=Cj->child+Cj->numChilds; ++Ci) {
      downwardPass(Ci);
    }
  }
}
#endif
