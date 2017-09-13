#ifndef helmholtz_h
#define helmholtz_h
#include <complex>
#include "exafmm.h"

namespace exafmm {
  void P2P(Cell * Ci, Cell * Cj) {
    Body * Bi = Ci->body;
    Body * Bj = Cj->body;
    for (int i=0; i<Ci->numBodies; i++) {
      complex_t p = 0;
      cvec3 F = complex_t(0.,0.);
      for (int j=0; j<Cj->numBodies; j++) {
        vec3 dX = Bi[i].X - Bj[j].X;
        real_t R2 = norm(dX);
        if (R2 != 0) {
          real_t R = std::sqrt(R2);
          complex_t pij = std::exp(I * R * WAVEK) * Bj[j].q / R;
          complex_t coef = (1/R2 - I*WAVEK/R) * pij;
          p += pij;
          for (int d=0; d<3; d++) {
            F[d] += coef * dX[d];
          }
        }
      }
      Bi[i].p += p;
      Bi[i].F += F;
    }
  }
}
#endif
