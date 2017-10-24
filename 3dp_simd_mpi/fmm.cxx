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
  Cells cells = buildTree(bodies);
  Bodies jbodies = bodies;
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
  horizontalPass(cells, jcells);
  stop("M2L & P2P");
  start("L2L & L2P");
  downwardPass(icells);
  stop("L2L & L2P");
  totalFMM += stop("Total FMM");
  Bodies bodies2;
  if ((IMAGES == 0) | (args.numBodies < 1000)) {
    start("Direct N-Body");
    const int numTargets = 10;
    jbodies = bodies;
    sampleBodies(bodies, numTargets);
    bodies2 = bodies;
    initTarget(bodies);
    for (int irank=0; irank<MPISIZE; irank++) {
      shiftBodies(jbodies);
      direct(bodies, jbodies);
    }
    stop("Direct N-Body");
  } else {
    start("Dipole correction");
    vec3 localDipole = 0, globalDipole = 0;
    for (size_t b=0; b<bodies.size(); b++) localDipole += bodies[b].X * bodies[b].q;
    MPI_Allreduce(&localDipole[0], &globalDipole[0], 3, MPI_REAL_T, MPI_SUM, MPI_COMM_WORLD);
    real_t coef = 4 * M_PI / (3 * CYCLE * CYCLE * CYCLE);
    int numBodies, numBodiesLocal = bodies.size();
    MPI_Allreduce(&numBodiesLocal, &numBodies, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    real_t dnorm = norm(globalDipole) / numBodies;
    for (size_t b=0; b<bodies.size(); b++) {
      bodies[b].p -= coef * dnorm / bodies[b].q;
      bodies[b].F -= globalDipole * coef;
    }
    stop("Dipole correction");

    print("Ewald Profiling");
    bodies2 = bodies;
    initTarget(bodies);
    jbodies = bodies;
    gatherBodies(jbodies);
    start("Build tree");
    Cells jcells2 = buildTree(jbodies);
    stop("Build tree");
    start("Wave part");
    wavePart(ibodies, jbodies);
    stop("Wave part");
    start("Real part");
    realPart(cells, jcells2);
    selfTerm(bodies);
    stop("Real part");
  }

  double pDif, pNrm;
  if (IMAGES == 0) {
    double pDifLocal = getDifScalar(bodies, bodies2);
    double pNrmLocal = getNrmScalar(bodies2);
    MPI_Allreduce(&pDifLocal, &pDif, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&pNrmLocal, &pNrm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  } else {
    double pSum, pSumLocal = getSumScalar(bodies);
    double pSum2, pSumLocal2 = getSumScalar(bodies2);
    MPI_Allreduce(&pSumLocal, &pSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&pSumLocal2, &pSum2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    pDif = (pSum - pSum2) * (pSum - pSum2);
    pNrm = pSum * pSum;
  }
  double pRel = std::sqrt(pDif/pNrm);
  double FDif, FDifLocal = getDifVector(bodies, bodies2);
  double FNrm, FNrmLocal = getNrmVector(bodies2);
  MPI_Allreduce(&FDifLocal, &FDif, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&FNrmLocal, &FNrm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  double FRel = std::sqrt(FDif/FNrm);
  print("FMM vs. direct");
  print("Rel. L2 Error (p)", pRel, false);
  print("Rel. L2 Error (F)", FRel, false);
  return 0;
}
