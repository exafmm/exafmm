#include "args.h"
#include "build_tree.h"
#include "dataset.h"
#include "local_essential_tree.h"
#include "partition.h"
#include "test.h"
using namespace exafmm;

int main(int argc, char ** argv) {
  Args args(argc, argv);
  NCRIT = args.ncrit;
  LEVEL = args.level;
  VERBOSE = args.verbose;

  Bodies bodies = initBodies(args.numBodies, args.distribution, MPIRANK, MPISIZE);
  for (size_t b=0; b<bodies.size(); b++) bodies[b].q = 1;

  partition(bodies);
  initKernel();
  Cells cells = buildTree(bodies);
  upwardPass(&cells[0]);
  localEssentialTree(bodies, cells);
  upwardPassLET(&cells[0]);

  int numBodies, numBodiesLocal = args.numBodies;
  MPI_Allreduce(&numBodiesLocal, &numBodies, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  print("numBodies", numBodies);
  print("cells[0].M[0]", cells[0].M[0]);
  assert(numBodies == std::real(cells[0].M[0]));
  print("Assertion passed");
  return 0;
}
