#ifndef traverse_eager_h
#define traverse_eager_h
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

  //! Recursive call to dual tree traversal for horizontal pass
  void horizontalPass(Cell * Ci, Cell * Cj) {
    vec2 dX = Ci->X - Cj->X;
    real_t R2 = norm(dX) * THETA * THETA;
    if (R2 > (Ci->R + Cj->R) * (Ci->R + Cj->R)) {
      M2L(Ci, Cj);
    } else if (Ci->numChilds == 0 && Cj->numChilds == 0) {
      P2P(Ci, Cj);
    } else if (Cj->numChilds == 0 || (Ci->R >= Cj->R && Ci->numChilds != 0)) {
      for (Cell * ci=Ci->child; ci!=Ci->child+Ci->numChilds; ci++) {
#pragma omp task untied if(ci->numBodies > 100)
        horizontalPass(ci, Cj);
      }
    } else {
      for (Cell * cj=Cj->child; cj!=Cj->child+Cj->numChilds; cj++) {
        horizontalPass(Ci, cj);
      }
    }
#pragma omp taskwait
  }

  //! Horizontal pass interface
  void horizontalPass(Cells & icells, Cells & jcells) {
#pragma omp parallel
#pragma omp single nowait
    horizontalPass(&icells[0], &jcells[0]);
  }

  //! Recursive call to pre-order tree traversal for downward pass
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
