#include "args.h"
#include "base_mpi.h"
#include "build_tree.h"
#include "dataset.h"
#include "kernel.h"
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

  Bodies bodies = initBodies(numBodies, distribution, baseMPI.mpirank, baseMPI.mpisize); // Initialize bodies
  real_t r0, x0[2], R0, X0[2];                                  // Initialize local & global bounds
  getBounds(bodies, r0, x0);                                    // Get local bounds
  allreduceBounds(R0, X0, r0, x0);                              // Reduce to global bounds
#if 0
  if (baseMPI.mpirank == 0) std::cout << "global R: " << R0 << std::endl;
#endif
  return 0;
}
