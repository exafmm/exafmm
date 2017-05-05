#define EXAFMM_LAZY 1
#include <cassert>
#include "build_tree.h"
#include "dataset.h"
#include "test.h"
#include "traverse_lazy.h"
using namespace exafmm;

int main(int argc, char ** argv) {
  const int numBodies = atoi(argv[1]);                          // Number of bodies
  theta = atof(argv[2]);                                        // Multipole acceptance criterion
  ncrit = atoi(argv[3]);                                        // Number of bodies per leaf cell
  const char * distribution = argv[4];                          // Type of distribution

  Bodies bodies = initBodies(numBodies, distribution);
  for (size_t b=0; b<bodies.size(); b++) {                      // Loop over bodies
    bodies[b].q = 1;                                            // Initialize with unit charge
  }

  Cells cells = buildTree(bodies);                              // Build tree

  test::upwardPass(&cells[0]);                                  // Upward pass

  getList(&cells[0], &cells[0]);                                // Create interaction list
  test::evaluate(cells);                                        // Horizontal pass

  test::downwardPass(&cells[0]);                                // Downward pass
  
  // Check answer
  printf("%-20s : %i\n", "Num of Bodies", numBodies);
  printf("--- %-18s ------------\n", "Checking potential");
  for (size_t b=0; b<bodies.size(); b++) {
    assert(numBodies == bodies[b].p);
  }
  printf("--- %-18s ------------\n", "Assertion Passed!");
  return 0;
}
