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
    Ci->M.resize(P, 0.0);
    Ci->L.resize(P, 0.0);
    if (Ci->numChilds == 0) P2M(Ci);
    M2M(Ci);
  }

  //! Upward pass interface
  void upwardPass(Cells & cells) {
#pragma omp parallel
#pragma omp single nowait
    upwardPass(&cells[0]);
  }

  //! Recursive call to dual tree traversal for list construction
  void getList(Cell * Ci, Cell * Cj) {
    vec2 dX = Ci->X - Cj->X;
    real_t R2 = norm(dX) * THETA * THETA;
    if (R2 > (Ci->R + Cj->R) * (Ci->R + Cj->R)) {
      Ci->listM2L.push_back(Cj);
    } else if (Ci->numChilds == 0 && Cj->numChilds == 0) {
      Ci->listP2P.push_back(Cj);
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
#pragma omp parallel for schedule(dynamic)
    for (size_t i=0; i<cells.size(); i++) {
      for (size_t j=0; j<cells[i].listM2L.size(); j++) {
        M2L(&cells[i],cells[i].listM2L[j]);
      }
      for (size_t j=0; j<cells[i].listP2P.size(); j++) {
        P2P(&cells[i],cells[i].listP2P[j]);
      }
    }
  }

  //! Horizontal pass interface
  void horizontalPass(Cells & icells, Cells & jcells) {
    getList(&icells[0], &jcells[0]);
    evaluate(icells);
  }

  //! Recursive call to pre-order traversal for downward pass
  void downwardPass(Cell * Cj) {
    L2L(Cj);
    if (Cj->numChilds == 0) L2P(Cj);
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
    P2P(Ci, Cj);
  }
}

#endif
