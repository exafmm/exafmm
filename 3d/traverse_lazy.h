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

  //! Recursive call to dual tree traversal for list construction
  void getList(Cell * Ci, Cell * Cj) {
    vec3 dX = Ci->X - Cj->X;
    real_t R2 = norm(dX) * theta * theta;
    if (R2 > (Ci->R + Cj->R) * (Ci->R + Cj->R)) {
      Ci->listM2L.push_back(Cj);
    } else if (Ci->NCHILD == 0 && Cj->NCHILD == 0) {
      Ci->listP2P.push_back(Cj);
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
    P2P(Ci, Cj);
  }
}
#endif
