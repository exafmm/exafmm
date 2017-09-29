#include "args.h"
#include "build_tree.h"
#include "dataset.h"
#include "test.h"
using namespace exafmm;

int main(int argc, char ** argv) {
  Args args(argc, argv);
  THETA = args.theta;
  NCRIT = args.ncrit;
  VERBOSE = args.verbose;
  const int numBodies = args.numBodies;
  const char * distribution = args.distribution;
  CYCLE = 2 * M_PI;
  IMAGES = args.images;

  Bodies bodies = initBodies(numBodies, distribution);
  for (size_t b=0; b<bodies.size(); b++) bodies[b].q = 1;

  Cells cells = buildTree(bodies);
  initKernel();
  upwardPass(cells);
  horizontalPass(cells, cells);
  downwardPass(cells);

  uint64_t imageBodies = std::pow(3,3*IMAGES) * bodies.size();
  print("numBodies", imageBodies);
  print("bodies[0].p", bodies[0].p);
  for (size_t b=0; b<bodies.size(); b++) assert(imageBodies == bodies[b].p);
  print("Assertion passed");
  return 0;
}
