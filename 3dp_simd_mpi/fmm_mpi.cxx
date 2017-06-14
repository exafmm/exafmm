#include "args.h"
#include "base_mpi.h"
#include "build_tree.h"
#include "dataset.h"
#include "kernel.h"
#include "partition.h"
#include "timer.h"
#if EXAFMM_EAGER
#include "traverse_eager.h"
#elif EXAFMM_LAZY
#include "traverse_lazy.h"
#endif
#include "verify.h"
using namespace exafmm;

int main(int argc, char ** argv) {
  Args args(argc, argv);                                        // Argument parser
  P = args.P;                                                   // Order of expansions
  theta = args.theta;                                           // Multipole acceptance criterion
  ncrit = args.ncrit;                                           // Number of bodies per leaf cell
  const int numBodies = args.numBodies;                         // Number of bodies
  const char * distribution = args.distribution;                // Type of distribution
  BaseMPI baseMPI;                                              // Initialize MPI environment
  Partition partition;

  real_t r0, x0[3];                                             // Initialize local & global bounds
  Bodies bodies = initBodies(numBodies, distribution, baseMPI.mpirank, baseMPI.mpisize); // Initialize bodies
  Bounds localBounds = getBounds(bodies, r0, x0);               // Get local bounds
  Bounds globalBounds = allreduceBounds(localBounds);           // Reduce to global bounds
  partition.bisection(bodies, globalBounds);

  return 0;
}
