#ifndef kernel_h
#define kernel_h
#include "exafmm.h"

namespace exafmm {
  //! P2P kernel between cells Ci and Cj
  void P2P(Cell * Ci, Cell * Cj) {
    Body * Bi = Ci->body;
    Body * Bj = Cj->body;
    for (int i=0; i<Ci->numBodies; i++) {
      real_t p = 0;
      vec2 F = 0;
      for (int j=0; j<Cj->numBodies; j++) {
        vec2 dX = Bi[i].X - Bj[j].X;
        real_t R2 = norm(dX);
        if (R2 != 0) {
          real_t invR = 1 / sqrt(R2);
          real_t logR = Bj[j].q * log(invR);
          p += logR;
          F += dX * Bj[j].q / R2;
        }
      }
      Bi[i].p += p;
      Bi[i].F -= F;
    }
  }

  //! P2M kernel for cell C
  void P2M(Cell * C) {
    for (Body * B=C->body; B!=C->body+C->numBodies; B++) {
      vec2 dX = B->X - C->X;
      complex_t Z(dX[0],dX[1]), powZ(1.0, 0.0);
      C->M[0] += B->q;
      for (int n=1; n<P; n++) {
        powZ *= Z / real_t(n);
        C->M[n] += powZ * B->q;
      }
    }
  }

  //! M2M kernel for one parent cell Ci
  void M2M(Cell * Ci) {
    for (Cell * Cj=Ci->child; Cj!=Ci->child+Ci->numChilds; Cj++) {
      vec2 dX = Cj->X - Ci->X;
      for (int k=0; k<P; k++) {
        complex_t Z(dX[0],dX[1]), powZ(1.0, 0.0);
        Ci->M[k] += Cj->M[k];
        for (int n=1; n<=k; n++) {
          powZ *= Z / real_t(n);
          Ci->M[k] += Cj->M[k-n] * powZ;
        }
      }
    }
  }

  //! M2L kernel between cells Ci and Cj
  void M2L(Cell * Ci, Cell * Cj) {
    vec2 dX = Ci->X - Cj->X;
    complex_t Z(dX[0],dX[1]), powZn(1.0, 0.0), powZnk(1.0, 0.0), invZ(powZn/Z);
    Ci->L[0] += -Cj->M[0] * log(Z);
    Ci->L[0] += Cj->M[1] * invZ;
    powZn = invZ;
    for (int k=2; k<P; k++) {
      powZn *= real_t(k-1) * invZ;
      Ci->L[0] += Cj->M[k] * powZn;
    }
    Ci->L[1] += -Cj->M[0] * invZ;
    powZn = invZ;
    for (int k=1; k<P; k++) {
      powZn *= real_t(k) * invZ;
      Ci->L[1] += -Cj->M[k] * powZn;
    }
    real_t Cnk = -1;
    for (int n=2; n<P; n++) {
      Cnk *= -1;
      powZnk *= invZ;
      powZn = Cnk * powZnk;
      for (int k=0; k<P; k++) {
        powZn *= real_t(n+k-1) * invZ;
        Ci->L[n] += Cj->M[k] * powZn;
      }
      powZnk *= real_t(n-1);
    }
  }

  //! L2L kernel for one parent cell Cj
  void L2L(Cell * Cj) {
    for (Cell * Ci=Cj->child; Ci<Cj->child+Cj->numChilds; Ci++) {
      vec2 dX = Ci->X - Cj->X;
      complex_t Z(dX[0],dX[1]);
      for (int l=0; l<P; l++) {
        complex_t powZ(1.0, 0.0);
        Ci->L[l] += Cj->L[l];
        for (int k=1; k<P-l; k++) {
          powZ *= Z / real_t(k);
          Ci->L[l] += Cj->L[l+k] * powZ;
        }
      }
    }
  }

  //! L2P kernel for cell C
  void L2P(Cell * C) {
    for (Body * B=C->body; B!=C->body+C->numBodies; B++) {
      vec2 dX = B->X - C->X;
      complex_t Z(dX[0],dX[1]), powZ(1.0, 0.0);
      B->p += std::real(C->L[0]);
      B->F[0] += std::real(C->L[1]);
      B->F[1] -= std::imag(C->L[1]);
      for (int n=1; n<P; n++) {
        powZ *= Z / real_t(n);
        B->p += std::real(C->L[n] * powZ);
        if (n < P-1) {
          B->F[0] += std::real(C->L[n+1] * powZ);
          B->F[1] -= std::imag(C->L[n+1] * powZ);
        }
      }
    }
  }
}

#endif
