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
  Args args(argc, argv);
  P = args.P;
  theta = args.theta;
  ncrit = args.ncrit;
  const int numBodies = args.numBodies;
  const char * distribution = args.distribution;

  if (args.verbose) {
    printf("--- %-20s --------\n", "FMM Parameter");
    args.print(20);
  }

  int repeat = 10;
  Verify verify(args.path);
  verify.verbose = args.verbose;
  bool isAccuracy = true;
  double totalFMM = 0;
  for (int t=0; t<repeat; t++) {
    printf("--- %-20s --------\n", "FMM Profiling");
    printf("--- %-17s %2d --------\n", "Time average loop", t);
    start("Initialize bodies");
    Bodies bodies = initBodies(numBodies, distribution);
    stop("Initialize bodies");
    start("Total FMM");
    start("Build tree");
    Cells cells = buildTree(bodies);
    stop("Build tree");
    start("P2M & M2M");
    upwardPass(cells);
    stop("P2M & M2M");
    start("M2L & P2P");
    horizontalPass(cells, cells);
    stop("M2L & P2P");
    start("L2L & L2P");
    downwardPass(cells);
    stop("L2L & L2P");
    totalFMM += stop("Total FMM");

    if (isAccuracy) {
      start("Direct N-Body");
      const int numTargets = 10;
      Bodies jbodies = bodies;
      if (numBodies > numTargets) {
        int stride = bodies.size() / numTargets;
        for (int b=0; b<numTargets; b++) {
          bodies[b] = bodies[b*stride];
        }
        bodies.resize(numTargets);
      }
      Bodies bodies2 = bodies;
      for (size_t b=0; b<bodies.size(); b++) {
        bodies[b].p = 0;
        bodies[b].F = 0;
      }
      direct(bodies, jbodies);
      stop("Direct N-Body");

      double pDif = verify.getDifScalar(bodies, bodies2);
      double pNrm = verify.getNrmScalar(bodies2);
      double pRel = std::sqrt(pDif/pNrm);
      double FDif = verify.getDifVector(bodies, bodies2);
      double FNrm = verify.getNrmVector(bodies2);
      double FRel = std::sqrt(FDif/FNrm);
      printf("--- %-20s --------\n", "FMM vs. direct");
      printf("%-20s : %7.4e\n","Rel. L2 Error (p)", pRel);
      printf("%-20s : %7.4e\n","Rel. L2 Error (F)", FRel);
      isAccuracy = false;
      printf("--- %-20s --------\n", "Accuracy regression");
      bool pass = verify.accuracyRegression(args.getKey(), pRel, FRel);
      if (pass) printf("Passed accuracy regression\n");
      else abort();
      if (args.accuracy) std::exit(0);
    }
  }
  double averageFMM = totalFMM / repeat;
  bool pass = verify.timeRegression(args.getKey(), averageFMM);
  if (pass) printf("Passed time regression\n");
  else abort();
  return 0;
}
