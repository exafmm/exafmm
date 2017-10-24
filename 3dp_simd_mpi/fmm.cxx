#include "mpi_utils.h"
#include "args.h"
#include "build_tree.h"
#include "dataset.h"
#include "ewald.h"
#include "kernel.h"
#include "local_essential_tree.h"
#include "partition.h"
#include "timer.h"
#include "traverse.h"
#include "verify.h"
using namespace exafmm;

int main(int argc, char ** argv) {
  Args args(argc, argv);
  P = args.P;
  THETA = args.theta;
  NCRIT = args.ncrit;
  LEVEL = args.level;
  VERBOSE = args.verbose;
  IMAGES = args.images;
  CYCLE = 2 * M_PI;
  KSIZE = 11;
  ALPHA = KSIZE / CYCLE;
  SIGMA = .25 / M_PI;
  CUTOFF = CYCLE / 2;

  print("FMM Parameter");
  args.show();

  double totalFMM = 0;
  print("FMM Profiling");
  start("Initialize bodies");
  Bodies ibodies = initBodies(args.numBodies, args.distribution, MPIRANK, MPISIZE);
  stop("Initialize bodies");
  start("Total FMM");
  start("Partition");
  partition(ibodies);
  stop("Partition");
  start("Precalculation");
  initKernel();
  stop("Precalculation");
  start("Build tree");
  Cells icells = buildTree(ibodies);
  Bodies jbodies = ibodies;
  Cells jcells = buildTree(jbodies);
  stop("Build tree");
  start("P2M & M2M");
  upwardPass(jcells);
  stop("P2M & M2M");
  start("Local essential tree");
  localEssentialTree(jbodies, jcells);
  upwardPassLET(jcells);
  stop("Local essential tree");
  start("M2L & P2P");
  horizontalPass(icells, jcells);
  stop("M2L & P2P");
  start("L2L & L2P");
  downwardPass(icells);
  stop("L2L & L2P");
  totalFMM += stop("Total FMM");
  Bodies ibodies2;
  if ((IMAGES == 0) | (args.numBodies < 1000)) {
    start("Direct N-Body");
    const int numTargets = 10;
    jbodies = ibodies;
    sampleBodies(ibodies, numTargets);
    ibodies2 = ibodies;
    initTarget(ibodies);
    direct(ibodies, jbodies);
    stop("Direct N-Body");
  } else {
    start("Dipole correction");
    vec3 localDipole = 0, globalDipole = 0;
    for (size_t b=0; b<ibodies.size(); b++) localDipole += ibodies[b].X * ibodies[b].q;
    MPI_Allreduce(&localDipole[0], &globalDipole[0], 3, MPI_REAL_T, MPI_SUM, MPI_COMM_WORLD);
    real_t coef = 4 * M_PI / (3 * CYCLE * CYCLE * CYCLE);
    for (size_t b=0; b<ibodies.size(); b++) {
      real_t dnorm = norm(globalDipole);
      ibodies[b].p -= coef * dnorm / ibodies.size() / ibodies[b].q;
      ibodies[b].F -= globalDipole * coef;
    }
    stop("Dipole correction");

    print("Ewald Profiling");
    start("Build tree");
    ibodies2 = ibodies;
    initTarget(ibodies);
    jbodies = ibodies;
    jcells = buildTree(jbodies);
    stop("Build tree");
    start("Wave part");
    wavePart(ibodies, jbodies);
    stop("Wave part");
    start("Real part");
    realPart(icells, jcells);
    selfTerm(ibodies);
    stop("Real part");
  }

  double pDif, pNrm;
  if (IMAGES == 0) {
    pDif = getDifScalar(ibodies, ibodies2);
    pNrm = getNrmScalar(ibodies2);
  } else {
    double pSum = getSumScalar(ibodies);
    double pSum2 = getSumScalar(ibodies2);
    pDif = (pSum - pSum2) * (pSum - pSum2);
    pNrm = pSum * pSum;
  }
  double pRel = std::sqrt(pDif/pNrm);
  double FDif = getDifVector(ibodies, ibodies2);
  double FNrm = getNrmVector(ibodies2);
  double FRel = std::sqrt(FDif/FNrm);
  print("FMM vs. direct");
  print("Rel. L2 Error (p)", pRel, false);
  print("Rel. L2 Error (F)", FRel, false);
  return 0;
}
