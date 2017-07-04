#define EXAFMM_LAZY 1
#include <cassert>
#include "args.h"
#include "build_tree.h"
#include "dataset.h"
#include "test.h"
#include "traverse_lazy.h"
using namespace exafmm;

int main(int argc, char ** argv) {
  Args args(argc, argv);                                        // Argument parser
  theta = args.theta;                                           // Multipole acceptance criterion
  ncrit = args.ncrit;                                           // Number of bodies per leaf cell
  const int numBodies = args.numBodies;                         // Number of bodies
  const char * distribution = args.distribution;                // Type of distribution

  Bodies bodies = initBodies(numBodies, distribution);          // Initialize bodies
  for (size_t b=0; b<bodies.size(); b++) bodies[b].q = 1;       // Set unit charge

  Cells cells = buildTree(bodies);                              // Build tree
  test::upwardPass(&cells[0]);                                  // Upward pass
  getList(&cells[0], &cells[0]);                                // Create interaction list
  test::evaluate(cells);                                        // Horizontal pass
  test::downwardPass(&cells[0]);                                // Downward pass
  
  // Check answer
  printf("%-20s : %i\n", "Num of Bodies", numBodies);
  printf("--- %-18s ------------\n", "Checking potential");
  for (size_t b=0; b<bodies.size(); b++) assert(numBodies == bodies[b].p);
  printf("--- %-18s ------------\n", "Assertion Passed!");
  return 0;
}
