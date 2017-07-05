#include "args.h"
#include "build_tree.h"
#include "dataset.h"
#include "kernel.h"
#include "ewald.h"
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
  cycle = 2 * M_PI;
  images = 4;
  ksize = 11;
  alpha = ksize / cycle;
  sigma = .25 / M_PI;
  cutoff = cycle / 2;

  if (args.verbose) {
    printf("--- %-20s --------\n", "FMM Parameter");
    args.print(20);
  }

  int repeat = 10;
  Verify verify;
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
    initKernel();
    upwardPass(cells);
    stop("P2M & M2M");
    start("M2L & P2P");
    horizontalPass(cells, cells);
    stop("M2L & P2P");
    start("L2L & L2P");
    downwardPass(cells);
    stop("L2L & L2P");
    start("Dipole correction");
    vec3 dipole = 0;
    for (size_t b=0; b<bodies.size(); b++) dipole += bodies[b].X * bodies[b].q;
    real_t coef = 4 * M_PI / (3 * cycle * cycle * cycle);
    for (size_t b=0; b<bodies.size(); b++) {
      real_t dnorm = norm(dipole);
      bodies[b].p -= coef * dnorm / bodies.size() / bodies[b].q;
      bodies[b].F -= dipole * coef;
    }
    stop("Dipole correction");
    totalFMM += stop("Total FMM");

    if (isAccuracy) {
      printf("--- %-20s -------\n", "Ewald Profiling");
      start("Build tree");
      Bodies bodies2 = bodies;
      for (size_t b=0; b<bodies.size(); b++) {
        bodies[b].p = 0;
        bodies[b].F = 0;
      }
      Bodies jbodies = bodies;
      Cells  jcells = buildTree(jbodies);
      stop("Build tree");
      start("Wave part");
      wavePart(bodies, jbodies);
      stop("Wave part");
      start("Real part");
      realPart(&cells[0], &jcells[0]);
      selfTerm(bodies);
      stop("Real part");

      double pSum = verify.getSumScalar(bodies);
      double pSum2 = verify.getSumScalar(bodies2);
      double pDif = (pSum - pSum2) * (pSum - pSum2);
      double pNrm = pSum * pSum;
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
