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

  //! Recursive call to post-order tree traversal for evaluating monopoles
  void upwardPass(Cell * Ci) {
    for (Cell * Cj=Ci->CHILD; Cj!=Ci->CHILD+Ci->NCHILD; ++Cj) {   // Loop over child cells
      test::upwardPass(Cj);                                           //  Recursive call for child cell
    }                                                             // End loop over child cells
    Ci->M.resize(1, 0.0);                                         // Allocate and initialize multipole coefs
    Ci->L.resize(1, 0.0);                                         // Allocate and initialize local coefs
    if (Ci->NCHILD==0) test::P2M(Ci);                                          // If Ci is leaf
    test::M2M(Ci);                                                        // If Ci is not leaf
  }

  //! Recursive call to post-order tree traversal for evaluating monopoles
  void downwardPass(Cell * Cj) {
    test::L2L(Cj); 
    if (Cj->NCHILD==0) test::L2P(Cj);
    for (Cell * Ci=Cj->CHILD; Ci!=Cj->CHILD+Cj->NCHILD; ++Ci) {   // Loop over child cells
      test::downwardPass(Ci);                                           //  Recursive call for child cell
    }                                                             // End loop over child cells
  }
#if EXAFMM_LAZY
  void evaluate(Cells & cells) {
    for (size_t i=0; i<cells.size(); i++) {                     // Loop over cells
      for (size_t j=0; j<cells[i].listM2L.size(); j++) {        //  Loop over M2L list
        test::M2L(&cells[i], cells[i].listM2L[j]);
      }                                                         //  End loop over M2L list
      for (size_t j=0; j<cells[i].listP2P.size(); j++) {        //  Loop over P2P list
        test::P2P(&cells[i], cells[i].listP2P[j]);                     //   P2P kernel
      }                                                         //  End loop over P2P list
    }                                                           // End loop over cells
  }
#endif
}
#endif
