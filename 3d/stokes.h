#ifndef stokes_h
#define stokes_h
#include "exafmm.h"

namespace exafmm {
  void P2P(Cell * Ci, Cell * Cj) {
    Body * Bi = Ci->body;
    Body * Bj = Cj->body;
    for (int i=0; i<Ci->numBodies; i++) {
      vec3 p = 0;
      for (int j=0; j<Cj->numBodies; j++) {
        vec3 dX = Bi[i].X - Bj[j].X;
        real_t R2 = norm(dX);
        if (R2 != 0) {
          real_t fdX = Bj[j].q[0] * dX[0] + Bj[j].q[1] * dX[1] + Bj[j].q[2] * dX[2];
          real_t invR = 1.0 / std::sqrt(R2);
          real_t invR3 = invR / R2;
          p += Bj[j].q * invR;
          p += dX * invR3 * fdX;
        }
      }
      Bi[i].p += p;
    }
  }
}
#endif
