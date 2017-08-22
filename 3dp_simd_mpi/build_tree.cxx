#include <cassert>
#include "args.h"
#include "build_tree.h"
#include "dataset.h"
#include "test.h"
using namespace exafmm;

int main(int argc, char ** argv) {
  Args args(argc, argv);
  NCRIT = args.ncrit;
  VERBOSE = args.verbose;
  const int numBodies = args.numBodies;
  const char * distribution = args.distribution;

  Bodies bodies = initBodies(numBodies, distribution);
  for (size_t b=0; b<bodies.size(); b++) bodies[b].q = 1;

  Cells cells = buildTree(bodies);
  upwardPass(&cells[0]);

  print("numBodies", numBodies);
  print("cells[0].M[0]", cells[0].M[0]);
  assert(numBodies == std::real(cells[0].M[0]));
  print("Assertion passed");
  return 0;
}
