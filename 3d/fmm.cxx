#include "build_tree.h"
#include "kernel.h"
#include "timer.h"
#include "traversal.h"
using namespace exafmm;

int main(int argc, char ** argv) {
  const int numBodies = 10000;                                  // Number of bodies
  P = 10;                                                       // Order of expansions
  ncrit = 64;                                                   // Number of bodies per leaf cell
  theta = 0.4;                                                  // Multipole acceptance criterion

  printf("--- %-16s ------------\n", "FMM Profiling");          // Start profiling
  //! Initialize bodies
  start("Initialize bodies");                                   // Start timer
  Bodies bodies(numBodies);                                     // Initialize bodies
  real_t average = 0;                                           // Average charge
  srand48(0);                                                   // Set seed for random number generator
  for (int b=0; b<int(bodies.size()); b++) {                    // Loop over bodies
    for (int d=0; d<3; d++) {                                   //  Loop over dimension
      bodies[b].X[d] = drand48() * 2 * M_PI - M_PI;             //   Initialize positions
    }                                                           //  End loop over dimension
    bodies[b].q = drand48() - .5;                               //  Initialize charge
    average += bodies[b].q;                                     //  Accumulate charge
    bodies[b].p = 0;                                            //  Clear potential
    for (int d=0; d<3; d++) bodies[b].F[d] = 0;                 //  Clear force
  }                                                             // End loop over bodies
  average /= bodies.size();                                     // Average charge
  for (int b=0; b<int(bodies.size()); b++) {                    // Loop over bodies
    bodies[b].q -= average;                                     // Charge neutral
  }                                                             // End loop over bodies
  stop("Initialize bodies");                                    // Stop timer

  //! Build tree
  start("Build tree");                                          // Start timer
  Cells  cells = buildTree(bodies);                             // Build tree
  stop("Build tree");                                           // Stop timer

  //! FMM evaluation
#pragma omp parallel                                            // Start OpenMP
#pragma omp single nowait                                       // Start OpenMP single region with nowait
  {
    start("Upward pass");                                       // Start timer
    initKernel();                                               // Initialize kernel
    upwardPass(&cells[0]);                                      // Upward pass for P2M, M2M
    stop("Upward pass");                                        // Stop timer
    start("Traversal");                                         // Start timer
    traversal(&cells[0], &cells[0]);                            // Traversal for M2L, P2P
    stop("Traversal");                                          // Stop timer
    start("Downward pass");                                     // Start timer
    downwardPass(&cells[0]);                                    // Downward pass for L2L, L2P
    stop("Downward pass");                                      // Stop timer
  }

  //! Direct N-Body
  start("Direct N-Body");                                       // Start timer
  const int numTargets = 10;                                    // Number of targets for checking answer
  Bodies jbodies = bodies;                                      // Save bodies in jbodies
  int stride = bodies.size() / numTargets;                      // Stride of sampling
  for (int b=0; b<numTargets; b++) {                            // Loop over target samples
    bodies[b] = bodies[b*stride];                               //  Sample targets
  }                                                             // End loop over target samples
  bodies.resize(numTargets);                                    // Resize bodies
  Bodies bodies2 = bodies;                                      // Backup bodies
  for (int b=0; b<int(bodies.size()); b++) {                    // Loop over bodies
    bodies[b].p = 0;                                            //  Clear potential
    for (int d=0; d<3; d++) bodies[b].F[d] = 0;                 //  Clear force
  }                                                             // End loop over bodies
  direct(bodies, jbodies);                                      // Direct N-Body
  stop("Direct N-Body");                                        // Stop timer

  //! Verify result
  real_t pDif = 0, pNrm = 0, FDif = 0, FNrm = 0;
  for (int b=0; b<int(bodies.size()); b++) {                    // Loop over bodies & bodies2
    pDif += (bodies[b].p - bodies2[b].p) * (bodies[b].p - bodies2[b].p);// Difference of potential
    pNrm += bodies2[b].p * bodies2[b].p;                        //  Value of potential
    FDif += (bodies[b].F[0] - bodies2[b].F[0]) * (bodies[b].F[0] - bodies2[b].F[0]) +// Difference of force
      (bodies[b].F[1] - bodies2[b].F[1]) * (bodies[b].F[1] - bodies2[b].F[1]) +// Difference of force
      (bodies[b].F[2] - bodies2[b].F[2]) * (bodies[b].F[2] - bodies2[b].F[2]);// Difference of force
    FNrm += bodies[b].F[0] * bodies[b].F[0] + bodies[b].F[1] * bodies[b].F[1] +// Value of force
      bodies[b].F[2] * bodies[b].F[2];
  }                                                             // End loop over bodies & bodies2
  printf("--- %-16s ------------\n", "FMM vs. direct");         // Print message
  printf("%-20s : %8.5e s\n","Rel. L2 Error (p)", sqrt(pDif/pNrm));// Print potential error
  printf("%-20s : %8.5e s\n","Rel. L2 Error (F)", sqrt(FDif/FNrm));// Print force error
  return 0;
}
