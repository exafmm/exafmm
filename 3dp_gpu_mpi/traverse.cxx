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
  CYCLE = 2 * M_PI;
  IMAGES = args.images;

  Bodies bodies = initBodies(args.numBodies, args.distribution, MPIRANK, MPISIZE);
  for (size_t b=0; b<bodies.size(); b++) bodies[b].q = 1;

  partition(bodies);
  initKernel();
  Cells cells = buildTree(bodies);
  Bodies jbodies = bodies;
  Cells jcells = buildTree(jbodies);
  upwardPass(jcells);
  localEssentialTree(jbodies, jcells);
  upwardPassLET(jcells);
  horizontalPass(cells, jcells);
  downwardPass(cells);

  int numBodies, numBodiesLocal = bodies.size();
  MPI_Allreduce(&numBodiesLocal, &numBodies, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  uint64_t imageBodies = std::pow(3,3*IMAGES) * numBodies;
  print("numBodies", imageBodies);
  print("bodies[0].p", bodies[0].p);
  for (size_t b=0; b<bodies.size(); b++) assert(imageBodies == bodies[b].p);
  print("Assertion passed");
  return 0;
}
