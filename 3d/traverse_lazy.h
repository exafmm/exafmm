#ifndef traverse_lazy_h
#define traverse_lazy_h
#include "exafmm.h"

namespace exafmm {
  //! Recursive call to post-order tree traversal for upward pass
  void upwardPass(Cell * Ci) {
    for (Cell * Cj=Ci->CHILD; Cj!=Ci->CHILD+Ci->NCHILD; Cj++) { // Loop over child cells
#pragma omp task untied if(Cj->NBODY > 100)                     //  Start OpenMP task if large enough task
      upwardPass(Cj);                                           //  Recursive call for child cell
    }                                                           // End loop over child cells
#pragma omp taskwait                                            // Synchronize OpenMP tasks
    Ci->M.resize(NTERM, 0.0);                                   // Allocate and initialize multipole coefs
    Ci->L.resize(NTERM, 0.0);                                   // Allocate and initialize local coefs
    if(Ci->NCHILD==0) P2M(Ci);                                  // P2M kernel
    M2M(Ci);                                                    // M2M kernel
  }

  //! Upward pass interface
  void upwardPass(Cells & cells) {
#pragma omp parallel                                            // Start OpenMP
#pragma omp single nowait                                       // Start OpenMP single region with nowait
    upwardPass(&cells[0]);                                      // Pass root cell to recursive call
  }

  //! Recursive call to dual tree traversal for list construction
  void getList(Cell * Ci, Cell * Cj) {
    for (int d=0; d<3; d++) dX[d] = Ci->X[d] - Cj->X[d];        // Distance vector from source to target
    real_t R2 = norm(dX) * theta * theta;                       // Scalar distance squared
    if (R2 > (Ci->R + Cj->R) * (Ci->R + Cj->R)) {               // If distance is far enough
      Ci->listM2L.push_back(Cj);                                //  Add to M2L list
    } else if (Ci->NCHILD == 0 && Cj->NCHILD == 0) {            // Else if both cells are leafs
      Ci->listP2P.push_back(Cj);                                //  Add to P2P list
    } else if (Cj->NCHILD == 0 || (Ci->R >= Cj->R && Ci->NCHILD != 0)) {// If Cj is leaf or Ci is larger
      for (Cell * ci=Ci->CHILD; ci!=Ci->CHILD+Ci->NCHILD; ci++) {// Loop over Ci's children
        getList(ci, Cj);                                        //   Recursive call to target child cells
      }                                                         //  End loop over Ci's children
    } else {                                                    // Else if Ci is leaf or Cj is larger
      for (Cell * cj=Cj->CHILD; cj!=Cj->CHILD+Cj->NCHILD; cj++) {// Loop over Cj's children
        getList(Ci, cj);                                        //   Recursive call to source child cells
      }                                                         //  End loop over Cj's children
    }                                                           // End if for leafs and Ci Cj size
  }

    //! Evaluate M2L, P2P kernels
  void evaluate(Cells & cells) {
#pragma omp parallel for schedule(dynamic)
    for (size_t i=0; i<cells.size(); i++) {                     // Loop over cells
      for (size_t j=0; j<cells[i].listM2L.size(); j++) {        //  Loop over M2L list
        M2L(&cells[i],cells[i].listM2L[j]);                     //   M2L kernel
      }                                                         //  End loop over M2L list
      for (size_t j=0; j<cells[i].listP2P.size(); j++) {        //  Loop over P2P list
        P2P(&cells[i],cells[i].listP2P[j]);                     //   P2P kernel
      }                                                         //  End loop over P2P list
    }                                                           // End loop over cells
  }

  //! Horizontal pass interface
  void horizontalPass(Cells & icells, Cells & jcells) {
    getList(&icells[0], &jcells[0]);                            // Pass root cell to recursive call
    evaluate(icells);                                           // Evaluate M2L & P2P kernels
  }

  //! Recursive call to pre-order tree traversal for downward pass
  void downwardPass(Cell * Cj) {
    L2L(Cj);                                                    // L2L kernel
    if (Cj->NCHILD==0) L2P(Cj);                                 // L2P kernel
    for (Cell * Ci=Cj->CHILD; Ci!=Cj->CHILD+Cj->NCHILD; Ci++) { // Loop over child cells
#pragma omp task untied if(Ci->NBODY > 100)                     //  Start OpenMP task if large enough task
      downwardPass(Ci);                                         //  Recursive call for child cell
    }                                                           // End loop over chlid cells
#pragma omp taskwait                                            // Synchronize OpenMP tasks
  }

  //! Downward pass interface
  void downwardPass(Cells & cells) {
#pragma omp parallel                                            // Start OpenMP
#pragma omp single nowait                                       // Start OpenMP single region with nowait
    downwardPass(&cells[0]);                                    // Pass root cell to recursive call
  }

  //! Direct summation
  void direct(Bodies & bodies, Bodies & jbodies) {
    Cells cells(2);                                             // Define a pair of cells to pass to P2P kernel
    Cell * Ci = &cells[0];                                      // Allocate single target
    Cell * Cj = &cells[1];                                      // Allocate single source
    Ci->BODY = &bodies[0];                                      // Iterator of first target body
    Ci->NBODY = bodies.size();                                  // Number of target bodies
    Cj->BODY = &jbodies[0];                                     // Iterator of first source body
    Cj->NBODY = jbodies.size();                                 // Number of source bodies
    P2P(Ci, Cj);                                                // Evaluate P2P kenrel
  }
}
#endif
