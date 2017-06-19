#include "args.h"
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
  Args args(argc, argv);                              // Argument parser
  P = args.P;                                         // Order of expansions
  theta = args.theta;                                 // Multipole acceptance criterion
  ncrit = args.ncrit;                                 // Number of bodies per leaf cell
  const int numBodies = args.numBodies;               // Number of bodies
  const char * distribution = args.distribution;      // Type of distribution

  // Print arguments 
  if (args.verbose) {
    printf("--- %-16s ------------\n", "FMM Parameter");
    args.print(20);
  }

  // FMM
  int repeat = 10;                                    // Number of runs for time regression test
  Verify verify;                                      // Initialize verify object
  verify.verbose = args.verbose;                      // Set verbosity
  bool isAccuracy = true;                             // Flag for checking accuracy
  double totalFMM = 0;                                // Initialize total FMM time
  for (int t=0; t<repeat; t++) {                      // Loop over identical runs for time regression
    printf("--- %-16s ------------\n", "FMM Profiling");
    printf("--- %-17s %d ---------\n", "Time average loop", t);
    // Initialize bodies
    start("Initialize bodies");                       //  Start timer
    Bodies bodies = initBodies(numBodies, distribution); //  Initialize bodies
    stop("Initialize bodies");                        //  Stop timer

    start("Total FMM");                               //  Start FMM timer
    // Build tree
    start("Build tree");                              //  Start timer
    Cells cells = buildTree(bodies);                  //  Build tree
    stop("Build tree");                               //  Stop timer

    // FMM evaluation
    start("P2M & M2M");                               //  Start timer
    upwardPass(cells);                                //  Upward pass for P2M, M2M
    stop("P2M & M2M");                                //  Stop timer
    start("M2L & P2P");                               //  Start timer
    horizontalPass(cells, cells);                     //  Horizontal pass for M2L, P2P
    stop("M2L & P2P");                                //  Stop timer
    start("L2L & L2P");                               //  Start timer
    downwardPass(cells);                              //  Downward pass for L2L, L2P
    stop("L2L & L2P");                                //  Stop timer
    totalFMM += stop("Total FMM");                    //  Stop FMM timer

    if (isAccuracy) {                                 //  If perform accuracy regression test
      // Direct N-Body
      start("Direct N-Body");                         //   Start timer
      const int numTargets = 10;                      //   Number of targets for checking answer
      Bodies jbodies = bodies;                        //   Save bodies in jbodies
      if (numBodies > numTargets) {                   //   If bodies are more than sampled targets
        int stride = bodies.size() / numTargets;      //    Stride of sampling
        for (int b=0; b<numTargets; b++) {            //    Loop over target samples
          bodies[b] = bodies[b*stride];               //     Sample targets
        }                                             //    End loop over target samples
        bodies.resize(numTargets);                    //    Resize bodies
      }                                               //   End if
      Bodies bodies2 = bodies;                        //   Backup bodies
      for (size_t b=0; b<bodies.size(); b++) {        //   Loop over bodies
        bodies[b].p = 0;                              //    Clear potential
        for (int d=0; d<2; d++) bodies[b].F[d] = 0;   //    Clear force
      }                                               //   End loop over bodies
      direct(bodies, jbodies);                        //   Direct N-Body
      stop("Direct N-Body");                          //   Stop timer

      // Verify result
      double pDif = verify.getDifScalar(bodies, bodies2);
      double pNrm = verify.getNrmScalar(bodies2);
      double pRel = std::sqrt(pDif/pNrm);
      double FDif = verify.getDifVector(bodies, bodies2);
      double FNrm = verify.getNrmVector(bodies2);
      double FRel = std::sqrt(FDif/FNrm);
      printf("--- %-16s ------------\n", "FMM vs. direct");
      printf("%-20s : %8.5e s\n","Rel. L2 Error (p)", pRel); //   Print potential error
      printf("%-20s : %8.5e s\n","Rel. L2 Error (F)", FRel); //   Print force error
      isAccuracy = false;                                    //   Set accuracy check flag to false
      printf("--- %-19s ---------\n", "Accuracy regression");
      bool pass = verify.accuracyRegression(args.getKey(), pRel, FRel); //   Accuracy regression test
      if (!pass) abort();
      if (args.accuracy) std::exit(0);                //   Finish execution if only needs accuracy regression
    }                                                 //  End if perform accuracy regression test
  }                                                   // End loop over identical runs for time regression
  double averageFMM = totalFMM / repeat;              // Average FMM time
  bool pass = verify.timeRegression(args.getKey(), averageFMM);         // Time regression test
  if (!pass) abort();
  return 0;
}
