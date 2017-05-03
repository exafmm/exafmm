#define EXAFMM_LAZY 1
#include <cassert>
#include "build_tree.h"
#include "test.h"
#include "traverse_lazy.h"
using namespace exafmm;

int main(int argc, char ** argv) {
  const int numBodies = atoi(argv[1]);                          // Number of bodies
  ncrit = 64;                                                   // Number of bodies per leaf cell
  theta = 0.4;                                                  // Multipole acceptance criterion

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
  Cells cells = buildTree(bodies);

  //! Upward pass
  test::upwardPass(&cells[0]);

  //! Horizontal pass
  getList(&cells[0], &cells[0]);
  test::evaluate(cells);

  //! Downward pass
  test::downwardPass(&cells[0]);
  
  //! Check answer
  printf("%-20s : %i\n", "Num of Bodies", numBodies);
  printf("--- %-18s ------------\n", "Checking potential");
  for (size_t b=0; b<bodies.size(); b++) {
    assert(numBodies == bodies[b].p);
  }
  printf("--- %-18s ------------\n", "Assertion Passed!");
  return 0;
}
