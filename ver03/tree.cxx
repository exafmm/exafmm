#include <cassert>
#include "build_tree.h"
#include "test.h"
using namespace exafmm;

int main(int argc, char ** argv) {
  ncrit = 64;                                                   // Number of bodies per leaf cell
  const int numBodies = 100000;                                 // Number of bodies
  Bodies bodies(numBodies);
  for (size_t b=0; b<numBodies; ++b) {                          // Loop over bodies
    for (int d=0; d<2; d++) {                                   //  Loop over dimension
      bodies[b].X[d] = drand48();                               //   Initialize coordinates
      bodies[b].F[d] = 0;                                       //   Initialize force
    }                                                           //  End loop over dimension
    bodies[b].q = 1;                                            //  Initialize charges
    bodies[b].p = 0;                                            //  Initialize potential
  }                                                             // End loop over bodies

  Cells cells = buildTree(bodies);                              // Build tree

  test::upwardPass(&cells[0]);                                  // Upward pass

  // Check Answer
  printf("%-20s : %i\n", "Num of Bodies", numBodies);
  printf("--- %-18s ------------\n", "Checking monopole");
  assert(numBodies == std::real(cells[0].M[0]));
  printf("--- %-18s ------------\n", "Assertion Passed!");

  return 0;
}
