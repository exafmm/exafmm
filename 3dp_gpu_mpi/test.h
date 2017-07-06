#ifndef test_h
#define test_h
#include "exafmm.h"

namespace test {
  using exafmm::Cell;
  using exafmm::Body;
  using exafmm::Cells;
  using exafmm::Cells;

  void P2M(Cell * C) {
    for (Body * B=C->BODY; B!=C->BODY+C->NBODY; ++B) C->M[0] += B->q;
  }

  void M2M(Cell * Ci) {
    for (Cell * Cj=Ci->CHILD; Cj!=Ci->CHILD+Ci->NCHILD; ++Cj) Ci->M[0] += Cj->M[0];
  }

  inline void M2L(Cell * Ci, Cell * Cj) {
    Ci->L[0] += Cj->M[0];
  }

  void L2L(Cell * Cj) {
    for (Cell * Ci=Cj->CHILD; Ci!=Cj->CHILD+Cj->NCHILD; ++Ci) Ci->L[0] += Cj->L[0];
  }

  void L2P(Cell * C) {
    for (Body * B=C->BODY; B!=C->BODY+C->NBODY; ++B) B->p += std::real(C->L[0]);
  }

  void P2P(Cell * Ci, Cell * Cj) {
    Body * Bi = Ci->BODY;
    Body * Bj = Cj->BODY;
    int ni = Ci->NBODY;
    int nj = Cj->NBODY;
    for (int i=0; i<ni; i++) {
      for (int j=0; j<nj; j++) {
        Bi[i].p += Bj[j].q;
      }
    }
  }

  void upwardPass(Cell * Ci) {
    for (Cell * Cj=Ci->CHILD; Cj!=Ci->CHILD+Ci->NCHILD; ++Cj) {
      test::upwardPass(Cj);
    }
    Ci->M.resize(1, 0.0);
    Ci->L.resize(1, 0.0);
    if (Ci->NCHILD==0) test::P2M(Ci);
    test::M2M(Ci);
  }

  void downwardPass(Cell * Cj) {
    test::L2L(Cj);
    if (Cj->NCHILD==0) test::L2P(Cj);
    for (Cell * Ci=Cj->CHILD; Ci!=Cj->CHILD+Cj->NCHILD; ++Ci) {
      test::downwardPass(Ci);
    }
  }
#if EXAFMM_LAZY
  void evaluate(Cells & cells) {
    for (size_t i=0; i<cells.size(); i++) {
      for (size_t j=0; j<cells[i].listM2L.size(); j++) {
        test::M2L(&cells[i], cells[i].listM2L[j]);
      }
      for (size_t j=0; j<cells[i].listP2P.size(); j++) {
        test::P2P(&cells[i], cells[i].listP2P[j]);
      }
    }
  }
#endif
}
#endif
