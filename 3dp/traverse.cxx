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
  CYCLE = 2 * M_PI;
  IMAGES = args.images;

  Bodies bodies = initBodies(args.numBodies, args.distribution);
  for (size_t b=0; b<bodies.size(); b++) bodies[b].q = 1;

  initKernel();
  Cells cells = buildTree(bodies);
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
