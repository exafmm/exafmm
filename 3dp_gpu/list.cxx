#define EXAFMM_LAZY 1
#include <cassert>
#include "args.h"
#include "build_tree.h"
#include "dataset.h"
#include "test.h"
#include "traverse_lazy.h"
using namespace exafmm;

int main(int argc, char ** argv) {
  Args args(argc, argv);
  theta = args.theta;
  ncrit = args.ncrit;
  const int numBodies = args.numBodies;
  const char * distribution = args.distribution;

  Bodies bodies = initBodies(numBodies, distribution);
  for (size_t b=0; b<bodies.size(); b++) bodies[b].q = 1;

  Cells cells = buildTree(bodies);
  test::upwardPass(&cells[0]);
  getList(&cells[0], &cells[0]);
  test::evaluate(cells);
  test::downwardPass(&cells[0]);

  printf("%-20s : %i\n", "Num of Bodies", numBodies);
  printf("--- %-18s ------------\n", "Checking potential");
  for (size_t b=0; b<bodies.size(); b++) assert(numBodies == bodies[b].p);
  printf("--- %-18s ------------\n", "Assertion Passed!");
  return 0;
}
