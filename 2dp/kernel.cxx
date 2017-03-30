#include "kernel.h"
using namespace exafmm;

int main(int argc, char ** argv) {
  P = atoi(argv[1]);

  // P2M
  Bodies jbodies(1);
  for (int d=0; d<2; d++) jbodies[0].X[d] = 2;
  jbodies[0].q = 1;
  Cells cells(4);
  Cell * Cj = &cells[0];
  Cj->X[0] = 3;
  Cj->X[1] = 1;
  Cj->R = 1;
  Cj->BODY = &jbodies[0];
  Cj->NBODY = jbodies.size();
  Cj->M.resize(P, 0.0);
  P2M(Cj);

  // M2M
  Cell * CJ = &cells[1];
  CJ->CHILD = Cj;
  CJ->NCHILD = 1;
  CJ->X[0] = 4;
  CJ->X[1] = 0;
  CJ->R = 2;
  CJ->M.resize(P, 0.0);
  M2M(CJ);

  // M2L
  Cell * CI = &cells[2];
  CI->X[0] = -4;
  CI->X[1] = 0;
  CI->R = 2;
  CI->L.resize(P, 0.0);
  M2L(CI, CJ);

  // L2L
  Cell * Ci = &cells[3];
  CI->CHILD = Ci;
  CI->NCHILD = 1;
  Ci->X[0] = -3;
  Ci->X[1] = 1;
  Ci->R = 1;
  Ci->L.resize(P, 0.0);
  L2L(CI);

  // L2P
  Bodies bodies(1);
  bodies[0].X[0] = -2;
  bodies[0].X[1] = 2;
  bodies[0].q = 1;
  bodies[0].p = 0;
  for (int d=0; d<2; d++) bodies[0].F[d] = 0;
  Ci->BODY = &bodies[0];
  Ci->NBODY = bodies.size();
  L2P(Ci);

  // P2P
  Bodies bodies2(1);
  for (size_t b=0; b<bodies2.size(); b++) {
    bodies2[b] = bodies[b];
    bodies2[b].p = 0;
    for (int d=0; d<2; d++) bodies2[b].F[d] = 0;
  }
  Cj->NBODY = jbodies.size();
  Ci->NBODY = bodies2.size();
  Ci->BODY = &bodies2[0];
  P2P(Ci, Cj);

  // Verify results
  real_t pDif = 0, pNrm = 0, FDif = 0, FNrm = 0;
  for (size_t b=0; b<bodies.size(); b++) {
    pDif += (bodies[b].p - bodies2[b].p) * (bodies[b].p - bodies2[b].p);
    pNrm += bodies[b].p * bodies[b].p;
    FDif += (bodies[b].F[0] - bodies2[b].F[0]) * (bodies[b].F[0] - bodies2[b].F[0]) +
      (bodies[b].F[1] - bodies2[b].F[1]) * (bodies[b].F[1] - bodies2[b].F[1]);
    FNrm += bodies[b].F[0] * bodies[b].F[0] + bodies[b].F[1] * bodies[b].F[1];
  }
  printf("%-20s : %8.5e s\n","Rel. L2 Error (p)", sqrt(pDif/pNrm));
  printf("%-20s : %8.5e s\n","Rel. L2 Error (F)", sqrt(FDif/FNrm));
  return 0;
}
