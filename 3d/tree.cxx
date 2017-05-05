#include <cassert>
#include "build_tree.h"
#include "dataset.h"
#include "test.h"
using namespace exafmm;

int main(int argc, char ** argv) {
  const int numBodies = atoi(argv[1]);                          // Number of bodies
  ncrit = atoi(argv[2]);                                        // Number of bodies per leaf cell
  const char * distribution = argv[3];                          // Type of distribution

  Bodies bodies = initBodies(numBodies, distribution);          // Initialize bodies
  for (size_t b=0; b<bodies.size(); b++) {                      // Loop over bodies
    bodies[b].q = 1;                                            // Initialize with unit charge
  }

  Cells cells = buildTree(bodies);                              // Build tree

  test::upwardPass(&cells[0]);                                  // Upward pass
  
  // Check Answer
  printf("%-20s : %i\n", "Num of Bodies", numBodies);
  printf("--- %-18s ------------\n", "Checking monopole");
  assert(numBodies == std::real(cells[0].M[0]));
  printf("--- %-18s ------------\n", "Assertion Passed!");

  return 0;
}
