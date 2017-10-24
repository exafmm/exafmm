#include "args.h"
#include "build_tree.h"
#include "dataset.h"
#include "test.h"
using namespace exafmm;

int main(int argc, char ** argv) {
  Args args(argc, argv);
  NCRIT = args.ncrit;
  VERBOSE = args.verbose;

  Bodies bodies = initBodies(args.numBodies, args.distribution);
  for (size_t b=0; b<bodies.size(); b++) bodies[b].q = 1;

  initKernel();
  Cells cells = buildTree(bodies);
  upwardPass(cells);

  print("numBodies", bodies.size());
  print("cells[0].M[0]", cells[0].M[0]);
  assert(bodies.size() == std::real(cells[0].M[0]));
  print("Assertion passed");
  return 0;
}
