#ifndef helmholtz_h
#define helmholtz_h
#include "exafmm.h"

namespace exafmm {
  int NQUAD, NQUAD2;
  std::vector<real_t> XQUAD, XQUAD2;
  std::vector<real_t> WQUAD, WQUAD2;

  void cart2sph(const vec3 & dX, real_t & r, real_t & theta, real_t & phi) {
    r = sqrt(norm(dX));
    theta = r == 0 ? 0 : acos(dX[2] / r);
    phi = atan2(dX[1], dX[0]);
  }

  void sph2cart(real_t r, real_t theta, real_t phi, const vec3 & spherical, vec3 & cartesian) {
    real_t invSinTheta = theta == 0 ? 0 : 1 / std::sin(theta);
    real_t invR = r == 0 ? 0 : 1 / r;
    cartesian[0] = std::sin(theta) * std::cos(phi) * spherical[0]
      + std::cos(theta) * std::cos(phi) * invR * spherical[1]
      - std::sin(phi) * invR * invSinTheta * spherical[2];
    cartesian[1] = std::sin(theta) * std::sin(phi) * spherical[0]
      + std::cos(theta) * std::sin(phi) * invR * spherical[1]
      + std::cos(phi) * invR * invSinTheta * spherical[2];
    cartesian[2] = std::cos(theta) * spherical[0]
      - std::sin(theta) * invR * spherical[1];
  }

  void polynomial(real_t x, int n, real_t & pol, real_t & der, real_t & sum) {
    sum = 0.5 + x * x * 1.5;
    real_t pk = 1;
    real_t pkp1 = x;
    if (n < 2) {
      der = 0;
      sum = 0.5;
      if (n == 0) return;
      der = 1;
      sum += x * x * 1.5;
      return;
    }
    for (int k=1; k<n; k++) {
      real_t pkm1 = pk;
      pk = pkp1;
      pkp1 = ((2 * k + 1) * x * pk - k * pkm1) / (k + 1);
      sum += pkp1 * pkp1 * (k + 1.5);
    }
    pol = pkp1;
    der = n * (x * pkp1 - pk) / (x * x - 1);
  }

  void legendre(int nq, std::vector<real_t> & xq, std::vector<real_t> & wq) {
    real_t pol = 0, der, sum;
    real_t h = M_PI / (2 * nq);
    for (int i=1; i<=nq; i++) {
      xq[nq-i] = std::cos((2 * i - 1) * h);
    }
    xq[nq/2] = 0;
    for (int i=0; i<nq/2; i++) {
      real_t xk = xq[i];
      int ifout = 0;
      for (int k=0; k<10; k++) {
        polynomial(xk,nq,pol,der,sum);
        real_t delta = -pol / der;
        xk += delta;
        if (fabs(delta) < EPS) ifout++;
        if (ifout == 3) break;
      }
      xq[i] = xk;
      xq[nq-i-1] = -xk;
    }
    for (int i=0; i<(nq+1)/2; i++) {
      polynomial(xq[i],nq,pol,der,sum);
      wq[i] = 1 / sum;
      wq[nq-i-1] = wq[i];
    }
  }

  void initKernel() {
    XQUAD.resize(P);
    XQUAD2.resize(2*P);
    WQUAD.resize(P);
    WQUAD2.resize(2*P);
    NQUAD = fmax(6, P);
    NQUAD2 = fmax(6, 2*P);
    legendre(NQUAD, XQUAD, WQUAD);
    legendre(NQUAD2, XQUAD2, WQUAD2);
  }

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
