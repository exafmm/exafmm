#include "build_tree.h"
using namespace exafmm;

//! Recursive call to post-order tree traversal for evaluating monopoles
void evalMonopole(Cell * Ci) {
  for (Cell * Cj=Ci->CHILD; Cj!=Ci->CHILD+Ci->NCHILD; ++Cj) {   // Loop over child cells
    evalMonopole(Cj);                                           //  Recursive call for child cell
  }                                                             // End loop over child cells
  Ci->M.resize(1, 0.0);                                         // Allocate and initialize multipole coefs
  if (Ci->NCHILD==0) {                                          // If Ci is leaf
    for (Body * B=Ci->BODY; B!=Ci->BODY+Ci->NBODY; ++B)         //  Loop over bodies in Ci
      Ci->M[0] += B->q;                                         //   Calculate monopole
  }                                                             // End if
  else {                                                        // If Ci is not leaf
    for (Cell * Cj=Ci->CHILD; Cj!=Ci->CHILD+Ci->NCHILD; ++Cj)   //  Loop over Ci's child cells
      Ci->M[0] += Cj->M[0];                                     //   Calculate monopole
  }                                                             // End else
}

int main(int argc, char ** argv) {
  const int numBodies = atoi(argv[1]);                          // Number of bodies
  ncrit = 64;                                                   // Number of bodies per leaf cell

  //! Initialize bodies
  Bodies bodies(numBodies);                                     // Initialize bodies
  srand48(0);                                                   // Set seed for random number generator
  for (size_t b=0; b<bodies.size(); b++) {                      // Loop over bodies
    for (int d=0; d<3; d++) {                                   //  Loop over dimension
      bodies[b].X[d] = drand48() * 2 * M_PI - M_PI;             //   Initialize positions
    }                                                           //  End loop over dimension
    bodies[b].q = 1;                                            //  Initialize with unit charge
    bodies[b].p = 0;                                            //  Clear potential
    for (int d=0; d<3; d++) bodies[b].F[d] = 0;                 //  Clear force
  }                                                             // End loop over bodies

  //! Build tree
  Cells cells = buildTree(bodies);                              // Build tree

  //! Evaluate monopoles (Upward Pass)
  evalMonopole(&cells[0]);
  
  printf("%-20s : %i\n", "Num of Bodies", numBodies);
  printf("%-20s : %f \n", "Root's Monopole", cells[0].M[0]);
  return 0;
}
