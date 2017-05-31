#include "build_tree.h"
#include "kernel.h"
#include "timer.h"
#include "traverse.h"
#include "verify.h"
using namespace exafmm;

int main(int argc, char ** argv) {
  P = 10;                                                       // Order of expansions
  theta = .4;                                                   // Multipole acceptance criterion
  ncrit = 64;                                                   // Number of bodies per leaf cell
  const int numBodies = 10000;                                  // Number of bodies
  // FMM
  Verify verify;
  verify.verbose = 1;                                           // Initialize a verify instance
  double totalFMM = 0;                                          // Initialize total FMM time
  printf("--- %-16s ------------\n", "FMM Profiling");          // Start profiling
  // Initialize bodies
  start("Initialize bodies");
  Bodies bodies(numBodies);
  for (size_t b=0; b<numBodies; ++b) {                          // Loop over bodies
    for (int d=0; d<2; d++) {                                   //  Loop over dimension
      bodies[b].X[d] = drand48();                               //   Initialize coordinates
      bodies[b].F[d] = 0;                                       //   Initialize force
    }                                                           //  End loop over dimension
    bodies[b].q = 1;                                            //  Initialize charges
    bodies[b].p = 0;                                            //  Initialize potential
  }                                                             // End loop over bodies
  stop("Initialize bodies");                                    // Stop timer

  start("Total FMM");
  // Build tree
  start("Build tree");                                          // Start timer
  Cells cells = buildTree(bodies);                              // Build tree
  stop("Build tree");                                           // Stop timer

  // FMM evaluation
  start("P2M & M2M");                                           // Start timer
  upwardPass(cells);                                            // Upward pass for P2M, M2M
  stop("P2M & M2M");                                            // Stop timer
  start("M2L & P2P");                                           // Start timer
  horizontalPass(cells, cells);                                 // Horizontal pass for M2L, P2P
  stop("M2L & P2P");                                            // Stop timer
  start("L2L & L2P");                                           // Start timer
  downwardPass(cells);                                          // Downward pass for L2L, L2P
  stop("L2L & L2P");                                            // Stop timer

  totalFMM += stop("Total FMM");

  // Direct N-Body
  start("Direct N-Body");                                       // Start timer
  const int numTargets = 10;                                    // Number of targets for checking answer
  Bodies jbodies = bodies;                                      // Save bodies in jbodies
  if (numBodies > numTargets) {                                 // If bodies are more than sampled targets
    int stride = bodies.size() / numTargets;                    //  Stride of sampling
    for (int b=0; b<numTargets; b++) {                          //  Loop over target samples
      bodies[b] = bodies[b*stride];                             //   Sample targets
    }                                                           //  End loop over target samples
    bodies.resize(numTargets);                                  //  Resize bodies
  }                                                             // End if
  Bodies bodies2 = bodies;                                      // Backup bodies
  for (size_t b=0; b<bodies.size(); b++) {                      // Loop over bodies
    bodies[b].p = 0;                                            //  Clear potential
    for (int d=0; d<2; d++) bodies[b].F[d] = 0;                 //  Clear force
  }                                                             // End loop over bodies
  direct(bodies, jbodies);                                      // Direct N-Body
  stop("Direct N-Body");                                        // Stop timer

  //! Verify result
  double pDif = verify.getDifScalar(bodies, bodies2);
  double pNrm = verify.getNrmScalar(bodies2);
  double pRel = std::sqrt(pDif/pNrm);
  double FDif = verify.getDifVector(bodies, bodies2);
  double FNrm = verify.getNrmVector(bodies2);
  double FRel = std::sqrt(FDif/FNrm);
  printf("--- %-16s ------------\n", "FMM vs. direct");          //  Print message
  printf("%-20s : %8.5e s\n","Rel. L2 Error (p)", pRel);         // Print potential error
  printf("%-20s : %8.5e s\n","Rel. L2 Error (F)", FRel);          // Print force error
  return 0;
}
