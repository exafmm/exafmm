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
  Coefs MJ, Mj, LI, Li;

  Kernel kernel(args.P, eps2, wavek);
  logger::verbose = true;

  Cells cells(4);
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



  C_iter Cj = cells.begin();
  Cj->X = 1;
  Cj->X[0] = 3;
  Cj->R = 1;
  Cj->S_BODY = jbodies.begin();
  Cj->S_NBODY = jbodies.size();
  Mj.resize(kernel.NTERM, 0.0);
  kernel.P2M(Cj, Mj);

#if 0
  C_iter CJ = cells.begin()+1;
  CJ->ICHILD = Cj-cells.begin();
  CJ->NCHILD = 1;
  CJ->X = 0;
  CJ->X[0] = 4;
  CJ->R = 2;
  MJ.resize(kernel.NTERM, 0.0);
  kernel.M2M(CJ, cells.begin(), MJ, Mj);

  C_iter CI = cells.begin()+2;
  CI->X = 0;
  CI->X[0] = -4;
  CI->R = 2;
  LI.resize(kernel.NTERM, 0.0);
  kernel.M2L(CI, CJ, MJ, LI);

  C_iter Ci = cells.begin()+3;
  Ci->X = 1;
  Ci->X[0] = -3;
  Ci->R = 1;
  Ci->IPARENT = 2;
  Li.resize(kernel.NTERM, 0.0);
  kernel.L2L(Ci, cells.begin(), Li, LI);
#else
  C_iter Ci = cells.begin()+3;
  Ci->X = 1;
  Ci->X[0] = -3;
  Ci->R = 1;
  Li.resize(kernel.NTERM, 0.0);
  kernel.M2L(Ci, Cj, Mj, Li);
#endif

  bodies[0].X = 2;
  bodies[0].X[0] = -2;
  //bodies[0].Q = 1;
  bodies[0].F = 0;
  Ci->T_BODY = bodies.begin();
  Ci->T_NBODY = bodies.size();
  kernel.L2P(Ci, Li);


  for (T_iter B=bodies2.begin(); B!=bodies2.end(); B++) {
    *B = bodies[B-bodies2.begin()];
    B->F = 0;
  }
  Cj->S_NBODY = jbodies.size();
  Ci->T_NBODY = bodies2.size();
  Ci->T_BODY = bodies2.begin();
  kernel.P2P(Ci, Cj);

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
