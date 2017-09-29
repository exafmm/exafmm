#include "args.h"
#include "build_tree.h"
#include "dataset.h"
#include "kernel.h"
#include "timer.h"
#include "traverse.h"
#include "verify.h"
using namespace exafmm;

int main(int argc, char ** argv) {
  Args args(argc, argv);
  P = args.P;
  THETA = args.theta;
  NCRIT = args.ncrit;
  VERBOSE = args.verbose;
  print("FMM Parameter");
  args.show();

  double totalFMM = 0;
  print("FMM Profiling");
  start("Initialize bodies");
  Bodies bodies = initBodies(args.numBodies, args.distribution);
  stop("Initialize bodies");
  start("Total FMM");
  start("Build tree");
  initKernel();
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

  start("Direct N-Body");
  const int numTargets = 10;
  Bodies jbodies = bodies;
  sampleBodies(bodies, numTargets);
  Bodies bodies2 = bodies;
  initTarget(bodies);
  direct(bodies, jbodies);
  stop("Direct N-Body");

  double pDif = getDifScalar(bodies, bodies2);
  double pNrm = getNrmScalar(bodies2);
  double pRel = std::sqrt(pDif/pNrm);
  double FDif = getDifVector(bodies, bodies2);
  double FNrm = getNrmVector(bodies2);
  double FRel = std::sqrt(FDif/FNrm);
  print("FMM vs. direct");
  print("Rel. L2 Error (p)", pRel, false);
  print("Rel. L2 Error (F)", FRel, false);
  return 0;
}
