#include "args.h"
#include "build_tree.h"
#include "dataset.h"
#include "local_essential_tree.h"
#include "partition.h"
#include "test.h"
using namespace exafmm;

int main(int argc, char ** argv) {
  Args args(argc, argv);
  THETA = args.theta;
  NCRIT = args.ncrit;
  LEVEL = args.level;
  VERBOSE = args.verbose;
  const int numBodies = args.numBodies;
  const char * distribution = args.distribution;
  CYCLE = 2 * M_PI;
  IMAGES = args.images;

  Bodies ibodies = initBodies(numBodies, distribution, MPIRANK, MPISIZE);
  for (size_t b=0; b<ibodies.size(); b++) ibodies[b].q = 1;

  partition(ibodies);
  initKernel();
  Cells icells = buildTree(ibodies);
  Bodies jbodies = ibodies;
  Cells jcells = buildTree(jbodies);
  upwardPass(jcells);
  localEssentialTree(jbodies, jcells);
  upwardPassLET(jcells);
  horizontalPass(icells, jcells);
  downwardPass(icells);

  uint64_t imageBodies = std::pow(3,3*IMAGES) * numBodies * MPISIZE;
  print("numBodies", imageBodies);
  print("bodies[0].p", ibodies[0].p);
  for (size_t b=0; b<ibodies.size(); b++) assert(imageBodies == ibodies[b].p);
  print("Assertion passed");
  return 0;
}
