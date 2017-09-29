#include "args.h"
#include "build_tree.h"
#include "dataset.h"
#include "test.h"
using namespace exafmm;

int main(int argc, char ** argv) {
  Args args(argc, argv);
  P = 1;
  THETA = args.theta;
  NCRIT = args.ncrit;
  VERBOSE = args.verbose;
  const int numBodies = args.numBodies;
  const char * distribution = args.distribution;

  Bodies bodies = initBodies(numBodies, distribution);
  for (size_t b=0; b<bodies.size(); b++) bodies[b].q = 1;

  Cells cells = buildTree(bodies);
  upwardPass(cells);
  horizontalPass(cells, cells);
  downwardPass(cells);

  print("numBodies", bodies.size());
  print("bodies[0].p", bodies[0].p);
  for (size_t b=0; b<bodies.size(); b++) assert(bodies.size() == bodies[b].p);
  print("Assertion passed");
  return 0;
}
