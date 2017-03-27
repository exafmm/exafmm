#ifndef traversal_h
#define traversal_h
#include "types.h"

namespace exafmm {
  int images;                                                   //!< Number of periodic image sublevels
  real_t theta;                                                 //!< Multipole acceptance criterion

  //! Recursive call for upward pass
  void upwardPass(Cell * Ci) {
    for (Cell * Cj=Ci->CHILD; Cj!=Ci->CHILD+Ci->NCHILD; Cj++) { // Loop over child cells
      upwardPass(Cj);                                           //  Recursive call
    }                                                           // End loop over child cells
    Ci->M.resize(P, 0);                                         // Allocate and initialize multipole coefs
    Ci->L.resize(P, 0);                                         // Allocate and initialize local coefs
    if (Ci->NCHILD == 0) P2M(Ci);                               // P2M kernel
    M2M(Ci);                                                    // M2M kernel
  }

  //! Dual tree traversal for a single pair of cells
  void traversal(Cell * Ci, Cell * Cj) {
    real_t dX[2];                                               // Distance vector
    for (int d=0; d<2; d++) dX[d] = Ci->X[d] - Cj->X[d] - Xperiodic[d];// Distance vector from source to target
    real_t R2 = norm(dX) * theta * theta;                       // Scalar distance squared
    if (R2 > (Ci->R + Cj->R) * (Ci->R + Cj->R)) {               // If distance is far enough
      M2L(Ci, Cj);                                              //  M2L kernel
    } else if (Ci->NCHILD == 0 && Cj->NCHILD == 0) {            // Else if both cells are leafs
      P2P(Ci, Cj);                                              //  P2P kernel
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

  //! Tree traversal of periodic cells
  void traversePeriodic(Cell * Ci0, Cell * Cj0, real_t cycle) {
    Cells pcells(9);                                            // Create cells
    for (int c=0; c<int(pcells.size()); c++) {                  // Loop over periodic cells
      pcells[c].M.resize(P, 0.0);                               //  Allocate & initialize M coefs
      pcells[c].L.resize(P, 0.0);                               //  Allocate & initialize L coefs
    }                                                           // End loop over periodic cells
    Cell * Ci = &pcells.back();                                 // Last cell is periodic parent cell
    *Ci = *Cj0;                                                 // Copy values from source root
    Ci->NCHILD = 8;                                             // Number of child cells for periodic center cell
    for (int level=0; level<images-1; level++) {                // Loop over sublevels of tree
      for (int ix=-1; ix<=1; ix++) {                            //  Loop over x periodic direction
        for (int iy=-1; iy<=1; iy++) {                          //   Loop over y periodic direction
          if (ix != 0 || iy != 0) {                             //    If periodic cell is not at center
            for (int cx=-1; cx<=1; cx++) {                      //     Loop over x periodic direction (child)
              for (int cy=-1; cy<=1; cy++) {                    //      Loop over y periodic direction (child)
                Xperiodic[0] = (ix * 3 + cx) * cycle;           //       Coordinate offset for x periodic direction
                Xperiodic[1] = (iy * 3 + cy) * cycle;           //       Coordinate offset for y periodic direction
                M2L(Ci0, Ci);                                   //       Perform M2L kernel
              }                                                 //      End loop over y periodic direction (child)
            }                                                   //     End loop over x periodic direction (child)
          }                                                     //    Endif for periodic center cell
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      Cell * Cj = &pcells.front();                              //  Iterator of periodic neighbor cells
      Ci->CHILD = Cj;                                           //  Child cells for periodic center cell
      for (int ix=-1; ix<=1; ix++) {                            //  Loop over x periodic direction
        for (int iy=-1; iy<=1; iy++) {                          //   Loop over y periodic direction
          if( ix != 0 || iy != 0) {                             //    If periodic cell is not at center
            Cj->X[0] = Ci->X[0] + ix * cycle;                   //     Set new x coordinate for periodic image
            Cj->X[1] = Ci->X[1] + iy * cycle;                   //     Set new y cooridnate for periodic image
            for (int n=0; n<P; n++) Cj->M[n] = Ci->M[n];        //     Copy multipoles to new periodic image
            Cj++;                                               //     Increment periodic cell iterator
          }                                                     //    Endif for periodic center cell
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      M2M(Ci);                                                  //  Evaluate periodic M2M kernels for this sublevel
      cycle *= 3;                                               //  Increase center cell size three times
    }                                                           // End loop over sublevels of tree
  }

  //! Evaluate P2P and M2L using dual tree traversal
  void traversal(Cell * Ci0, Cell * Cj0, real_t cycle) {
    if (images == 0) {                                          // If non-periodic boundary condition
      for (int d=0; d<2; d++) Xperiodic[d] = 0;                 //  No periodic shift
      traversal(Ci0, Cj0);                                      //  Traverse the tree
    } else {                                                    // If periodic boundary condition
      for (int ix=-1; ix<=1; ix++) {                            //  Loop over x periodic direction
        for (int iy=-1; iy<=1; iy++) {                          //   Loop over y periodic direction
          Xperiodic[0] = ix * cycle;                            //    Coordinate shift for x periodic direction
          Xperiodic[1] = iy * cycle;                            //    Coordinate shift for y periodic direction
          traversal(Ci0, Cj0);                                  //    Traverse the tree for this periodic image
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      traversePeriodic(Ci0, Cj0, cycle);                        //  Traverse tree for periodic images
    }                                                           // End if for periodic boundary condition
  }                                                             // End if for empty cell vectors

  //! Recursive call for downward pass
  void downwardPass(Cell * Cj) {
    L2L(Cj);                                                    // L2L kernel
    if (Cj->NCHILD == 0) L2P(Cj);                               // L2P kernel
    for (Cell * Ci=Cj->CHILD; Ci!=Cj->CHILD+Cj->NCHILD; Ci++) { // Loop over child cells
      downwardPass(Ci);                                         //  Recursive call
    }                                                           // End loop over child cells
  }

  //! Direct summation
  void direct(Bodies & bodies, Bodies & jbodies, real_t cycle) {
    Cells cells(2);                                             // Define a pair of cells to pass to P2P kernel
    Cell * Ci = &cells[0];                                      // Allocate single target cell
    Cell * Cj = &cells[1];                                      // Allocate single source cell
    Ci->BODY = &bodies[0];                                      // Pointer of first target body
    Ci->NBODY = bodies.size();                                  // Number of target bodies
    Cj->BODY = &jbodies[0];                                     // Pointer of first source body
    Cj->NBODY = jbodies.size();                                 // Number of source bodies
    int prange = 0;                                             // Range of periodic images
    for (int i=0; i<images; i++) {                              // Loop over periodic image sublevels
      prange += int(powf(3.,i));                                //  Accumulate range of periodic images
    }                                                           // End loop over perioidc image sublevels
    for (int ix=-prange; ix<=prange; ix++) {                    // Loop over x periodic direction
      for (int iy=-prange; iy<=prange; iy++) {                  //  Loop over y periodic direction
        Xperiodic[0] = ix * cycle;                              //   Coordinate shift for x periodic direction
        Xperiodic[1] = iy * cycle;                              //   Coordinate shift for y periodic direction
        P2P(Ci, Cj);                                            //   Evaluate P2P kernel
      }                                                         //  End loop over y periodic direction
    }                                                           // End loop over x periodic direction
  }
}

#endif
