#include "build_tree.h"
#include "kernel.h"
#include "ewald.h"
#include "timer.h"
#if EXAFMM_EAGER
#include "traverse_eager.h"
#elif EXAFMM_LAZY
#include "traverse_lazy.h"
#endif
using namespace exafmm;

int main(int argc, char ** argv) {
  const int numBodies = 1000;                                   // Number of bodies
  P = 10;                                                       // Order of expansions
  ncrit = 64;                                                   // Number of bodies per leaf cell
  cycle = 2 * M_PI;                                             // Cycle of periodic boundary condition
  theta = 0.4;                                                  // Multipole acceptance criterion
  images = 4;                                                   // 3^images * 3^images * 3^images periodic images

  ksize = 11;                                                   // Ewald wave number
  alpha = ksize / cycle;                                        // Ewald real/wave balance parameter
  sigma = .25 / M_PI;                                           // Ewald distribution parameter
  cutoff = cycle / 2;                                           // Ewald cutoff distance

  printf("--- %-16s ------------\n", "FMM Profiling");          // Start profiling
  //! Initialize bodies
  start("Initialize bodies");                                   // Start timer
  Bodies bodies(numBodies);                                     // Initialize bodies
  real_t average = 0;                                           // Average charge
  srand48(0);                                                   // Set seed for random number generator
  for (size_t b=0; b<bodies.size(); b++) {                      // Loop over bodies
    for (int d=0; d<3; d++) {                                   //  Loop over dimension
      bodies[b].X[d] = drand48() * cycle - cycle * .5;          //   Initialize positions
    }                                                           //  End loop over dimension
    bodies[b].q = drand48() - .5;                               //  Initialize charge
    average += bodies[b].q;                                     //  Accumulate charge
    bodies[b].p = 0;                                            //  Clear potential
    for (int d=0; d<3; d++) bodies[b].F[d] = 0;                 //  Clear force
  }                                                             // End loop over bodies
  average /= bodies.size();                                     // Average charge
  for (size_t b=0; b<bodies.size(); b++) {                      // Loop over bodies
    bodies[b].q -= average;                                     // Charge neutral
  }                                                             // End loop over bodies
  stop("Initialize bodies");                                    // Stop timer

  //! Build tree
  start("Build tree");                                          // Start timer
  Cells  cells = buildTree(bodies);                             // Build tree
  stop("Build tree");                                           // Stop timer

  //! FMM evaluation
  start("P2M & M2M");                                           // Start timer
  initKernel();                                                 // Initialize kernel
  upwardPass(cells);                                            // Upward pass for P2M, M2M
  stop("P2M & M2M");                                            // Stop timer
  start("M2L & P2P");                                           // Start timer
  horizontalPass(cells, cells);                                 // Horizontal pass for M2L, P2P
  stop("M2L & P2P");                                            // Stop timer
  start("L2L & L2P");                                           // Start timer
  downwardPass(cells);                                          // Downward pass for L2L, L2P
  stop("L2L & L2P");                                            // Stop timer

  //! Dipole correction
  start("Dipole correction");                                   // Start timer
  real_t dipole[3] = {0, 0, 0};                                 // Initialize dipole
  for (size_t b=0; b<bodies.size(); b++) {                      // Loop over bodies
    for (int d=0; d<3; d++) dipole[d] += bodies[b].X[d] * bodies[b].q;// Accumulate dipole
  }                                                             // End loop over bodies
  real_t coef = 4 * M_PI / (3 * cycle * cycle * cycle);         // Domain coefficient
  for (size_t b=0; b<bodies.size(); b++) {                      // Loop over bodies
    real_t dnorm = dipole[0] * dipole[0] + dipole[1] * dipole[1] + dipole[2] * dipole[2];// Norm of dipole
    bodies[b].p -= coef * dnorm / bodies.size() / bodies[b].q;  //  Correct potential
    for (int d=0; d!=3; d++) bodies[b].F[d] -= coef * dipole[d];//  Correct force
  }                                                             // End loop over bodies
  stop("Dipole correction");                                    // Stop timer

  printf("--- %-16s ------------\n", "Ewald Profiling");        // Print message
  //! Ewald summation
  start("Build tree");                                          // Start timer
  Bodies bodies2 = bodies;                                      // Backup bodies
  for (size_t b=0; b<bodies.size(); b++) {                      // Loop over bodies
    bodies[b].p = 0;                                            //  Clear potential
    for (int d=0; d<3; d++) bodies[b].F[d] = 0;                 //  Clear force
  }                                                             // End loop over bodies
  Bodies jbodies = bodies;                                      // Copy bodies
  Cells  jcells = buildTree(jbodies);                           // Build tree
  stop("Build tree");                                           // Stop timer
  start("Wave part");                                           // Start timer
  wavePart(bodies, jbodies);                                    // Ewald wave part
  stop("Wave part");                                            // Stop timer
  start("Real part");                                           // Start timer
  realPart(&cells[0], &jcells[0]);                              // Ewald real part
  selfTerm(bodies);                                             // Ewald self term
  stop("Real part");                                            // Stop timer

  //! Verify result
  real_t pSum = 0, pSum2 = 0, FDif = 0, FNrm = 0;
  for (size_t b=0; b<bodies.size(); b++) {                      // Loop over bodies & bodies2
    pSum += bodies[b].p * bodies[b].q;                          // Sum of potential for bodies
    pSum2 += bodies2[b].p * bodies2[b].q;                       // Sum of potential for bodies2
    FDif += (bodies[b].F[0] - bodies2[b].F[0]) * (bodies[b].F[0] - bodies2[b].F[0]) +// Difference of force
      (bodies[b].F[1] - bodies2[b].F[1]) * (bodies[b].F[1] - bodies2[b].F[1]) +// Difference of force
      (bodies[b].F[2] - bodies2[b].F[2]) * (bodies[b].F[2] - bodies2[b].F[2]);// Difference of force
    FNrm += bodies[b].F[0] * bodies[b].F[0] + bodies[b].F[1] * bodies[b].F[1] +// Value of force
      bodies[b].F[2] * bodies[b].F[2];
  }                                                             // End loop over bodies & bodies2
  real_t pDif = (pSum - pSum2) * (pSum - pSum2);                // Difference in sum
  real_t pNrm = pSum * pSum;                                    // Norm of the sum
  printf("--- %-16s ------------\n", "FMM vs. direct");         // Print message
  printf("%-20s : %8.5e s\n","Rel. L2 Error (p)", sqrt(pDif/pNrm));// Print potential error
  printf("%-20s : %8.5e s\n","Rel. L2 Error (F)", sqrt(FDif/FNrm));// Print force error
  return 0;
}
