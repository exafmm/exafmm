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
  
  Targets bodies(1), bodies2(1);
  Sources jbodies(1);
 
  vec3 XJ, Xj, XI, Xi;
  Coefs MJ, Mj, LI, Li;

  Kernel kernel(args.P, eps2, wavek);
  logger::verbose = true;

  Verify verify;
  jbodies[0].X = 2;
#if EXAFMM_BIOTSAVART
  jbodies[0].Q[0] = drand48();
  jbodies[0].Q[1] = drand48();
  jbodies[0].Q[2] = drand48();
  jbodies[0].Q[3] = 0.1;
#else
  jbodies[0].Q = 1;
#endif

  Xj = 1;
  Xj[0] = 3;
  Source* Bj = &jbodies[0];
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

  bodies[0].X = 2;
  bodies[0].X[0] = -2;
  bodies[0].F = 0;
 
  bodies2 = bodies;

  Target* Bi = &bodies[0];
  int ni = 1;
  kernel.L2P(Bi, ni, Xi, Li);

  Bi = &bodies2[0];
  kernel.P2P(Bi, ni, Bj, nj);

  std::fstream file;
  file.open("kernel.dat", std::ios::out | std::ios::app);
  double potDif = verify.getDifScalar(bodies, bodies2);
  double potNrm = verify.getNrmScalar(bodies);
  double accDif = verify.getDifVector(bodies, bodies2);
  double accNrm = verify.getNrmVector(bodies);
  std::cout << args.P << " " << std::sqrt(potDif/potNrm) << "  " << std::sqrt(accDif/accNrm) << std::endl;
  double potRel = std::sqrt(potDif/potNrm);
  double accRel = std::sqrt(accDif/accNrm);
  verify.print("Rel. L2 Error (pot)",potRel);
  verify.print("Rel. L2 Error (acc)",accRel);
  file << args.P << " " << std::sqrt(potDif/potNrm) << "  " << std::sqrt(accDif/accNrm) << std::endl;
  file.close();
  return 0;
}
