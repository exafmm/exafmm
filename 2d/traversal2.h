#ifndef traversal_h
#define traversal_h
#include "types.h"

namespace exafmm {
  real_t theta;                                                 //!< Multipole acceptance criterion
  std::multimap<Cell*,Cell*> listM2L;                           //!< M2L interaction list
  std::multimap<Cell*,Cell*> listP2P;                           //!< P2P interaction list
  typedef std::multimap<Cell*,Cell*>::iterator M_iter;          //!< Iterator of multimap

  //! Recursive call for upward pass
  void upwardPass(Cell * Ci) {
    for (Cell * Cj=Ci->CHILD; Cj!=Ci->CHILD+Ci->NCHILD; Cj++) { // Loop over child cells
#pragma omp task untied if(Cj->NBODY > 100)                     //  Start OpenMP task if large enough task
      upwardPass(Cj);                                           //  Recursive call
    }                                                           // End loop over child cells
#pragma omp taskwait                                            // Synchronize OpenMP tasks
    Ci->M.resize(P, 0);                                         // Allocate and initialize multipole coefs
    Ci->L.resize(P, 0);                                         // Allocate and initialize local coefs
    if (Ci->NCHILD == 0) P2M(Ci);                               // P2M kernel
    M2M(Ci);                                                    // M2M kernel
  }

  //! Dual tree traversal for a single pair of cells
  void traversal(Cell * Ci, Cell * Cj) {
    for (int d=0; d<2; d++) dX[d] = Ci->X[d] - Cj->X[d];        // Distance vector from source to target
    real_t R2 = norm(dX) * theta * theta;                       // Scalar distance squared
    if (R2 > (Ci->R + Cj->R) * (Ci->R + Cj->R)) {               // If distance is far enough
      listM2L.insert(std::pair<Cell*,Cell*>(Ci, Cj));           //  Add to M2L list
    } else if (Ci->NCHILD == 0 && Cj->NCHILD == 0) {            // Else if both cells are leafs
      listP2P.insert(std::pair<Cell*,Cell*>(Ci, Cj));           //  Add to P2P list
    } else if (Cj->NCHILD == 0 || (Ci->R >= Cj->R && Ci->NCHILD != 0)) {// Else if Cj is leaf or Ci is larger
      for (Cell * ci=Ci->CHILD; ci!=Ci->CHILD+Ci->NCHILD; ci++) {// Loop over Ci's children
        traversal(ci, Cj);                                      //   Traverse a single pair of cells
      }                                                         //  End loop over Ci's children
    } else {                                                    // Else if Ci is leaf or Cj is larger
      for (Cell * cj=Cj->CHILD; cj!=Cj->CHILD+Cj->NCHILD; cj++) {//  Loop over Cj's children
        traversal(Ci, cj);                                      //  Traverse a single pair of cells
      }                                                         //  End loop over Cj's children
    }                                                           // End if for leafs and Ci Cj size
  }

  //! Evaluate M2L, P2P kernels
  void evaluate() {
    for (M_iter M=listM2L.begin(); M!=listM2L.end(); M++) {     // Loop over all M2L lists
      M2L(M->first,M->second);                                  //  M2L kernel
    }                                                           // End loop over all M2L lists
    for (M_iter M=listP2P.begin(); M!=listP2P.end(); M++) {     // Loop over all P2P lists
      P2P(M->first,M->second);                                  //  P2P kernel
    }                                                           // End loop over all P2P lists
  }

  //! Recursive call for downward pass
  void downwardPass(Cell * Cj) {
    L2L(Cj);                                                    // L2L kernel
    if (Cj->NCHILD == 0) L2P(Cj);                               // L2P kernel
    for (Cell * Ci=Cj->CHILD; Ci!=Cj->CHILD+Cj->NCHILD; Ci++) { // Loop over child cells
#pragma omp task untied if(Ci->NBODY > 100)                     //  Start OpenMP task if large enough task
      downwardPass(Ci);                                         //  Recursive call
    }                                                           // End loop over child cells
#pragma omp taskwait                                            // Synchronize OpenMP tasks
  }

  //! Direct summation
  void direct(Bodies & bodies, Bodies & jbodies) {
    Cells cells(2);                                             // Define a pair of cells to pass to P2P kernel
    Cell * Ci = &cells[0];                                      // Allocate single target cell
    Cell * Cj = &cells[1];                                      // Allocate single source cell
    Ci->BODY = &bodies[0];                                      // Pointer of first target body
    Ci->NBODY = bodies.size();                                  // Number of target bodies
    Cj->BODY = &jbodies[0];                                     // Pointer of first source body
    Cj->NBODY = jbodies.size();                                 // Number of source bodies
    P2P(Ci, Cj);                                                // Evaluate P2P kernel
  }
}

#endif
