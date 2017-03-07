#include "args.h"
#include <fstream>
#include "kernel.h"
#include "namespace.h"
#include <vector>
#include "verify.h"
using namespace EXAFMM_NAMESPACE;

int main(int argc, char ** argv) {
  const real_t eps2 = 0.0;
  const complex_t wavek = complex_t(1.,.1) / real_t(2 * M_PI);
  Args args(argc, argv);
  
  Targets targets(1), targets2(1);
  Sources sources(1);
 
  vec3 XJ, Xj, XI, Xi;
  Coefs MJ, Mj, LI, Li;

  Kernel kernel(args.P, eps2, wavek);
  logger::verbose = true;

  Verify verify;
  sources[0].X = 2;
#if EXAFMM_BIOTSAVART
  sources[0].Q[0] = drand48();
  sources[0].Q[1] = drand48();
  sources[0].Q[2] = drand48();
  sources[0].Q[3] = 0.1;
#else
  sources[0].Q = 1;
#endif

  Xj = 1;
  Xj[0] = 3;
  Source* Bj = &sources[0];
  int nj = 1;
  Mj.resize(kernel.NTERM, 0.0);
  kernel.P2M(Xj, Mj, Bj, nj);

#if 1
  XJ = 0;
  XJ[0] = 4;
  MJ.resize(kernel.NTERM, 0.0);
  kernel.M2M(XJ, MJ, Xj, Mj);

  XI = 0;
  XI[0] = -4;
  LI.resize(kernel.NTERM, 0.0);
  kernel.M2L(XI, LI, XJ, MJ);

  Xi = 1;
  Xi[0] = -3;
  Li.resize(kernel.NTERM, 0.0);
  kernel.L2L(Xi, Li, XI, LI);
#else
  Xi = 1;
  Xi[0] = -3;
  Li.resize(kernel.NTERM, 0.0);
  kernel.M2L(Xi, Li, Xj, Mj);
#endif

  targets[0].X = 2;
  targets[0].X[0] = -2;
  targets[0].F = 0;
 
  targets2 = targets;

  Target* Bi = &targets[0];
  int ni = 1;
  kernel.L2P(Bi, ni, Xi, Li);

  Bi = &targets2[0];
  kernel.P2P(Bi, ni, Bj, nj);

  std::fstream file;
  file.open("kernel.dat", std::ios::out | std::ios::app);
  double potDif = verify.getDifScalar(targets, targets2);
  double potNrm = verify.getNrmScalar(targets);
  double accDif = verify.getDifVector(targets, targets2);
  double accNrm = verify.getNrmVector(targets);
  std::cout << args.P << " " << std::sqrt(potDif/potNrm) << "  " << std::sqrt(accDif/accNrm) << std::endl;
  double potRel = std::sqrt(potDif/potNrm);
  double accRel = std::sqrt(accDif/accNrm);
  verify.print("Rel. L2 Error (pot)",potRel);
  verify.print("Rel. L2 Error (acc)",accRel);
  file << args.P << " " << std::sqrt(potDif/potNrm) << "  " << std::sqrt(accDif/accNrm) << std::endl;
  file.close();
  return 0;
}
