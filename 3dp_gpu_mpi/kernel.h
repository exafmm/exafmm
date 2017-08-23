#ifndef kernel_h
#define kernel_h
#include "exafmm.h"

namespace exafmm {
  const complex_t I(0.,1.);

  inline int oddOrEven(int n) {
    return (((n) & 1) == 1) ? -1 : 1;
  }

  inline int ipow2n(int n) {
    return (n >= 0) ? 1 : oddOrEven(n);
  }

  void cart2sph(const vec3 & dX, real_t & r, real_t & theta, real_t & phi) {
    r = sqrt(norm(dX));
    theta = r == 0 ? 0 : acos(dX[2] / r);
    phi = atan2(dX[1], dX[0]);
  }

  void sph2cart(real_t r, real_t theta, real_t phi, const vec3 & spherical, vec3 & cartesian) {
    cartesian[0] = std::sin(theta) * std::cos(phi) * spherical[0]
      + std::cos(theta) * std::cos(phi) / r * spherical[1]
      - std::sin(phi) / r / std::sin(theta) * spherical[2];
    cartesian[1] = std::sin(theta) * std::sin(phi) * spherical[0]
      + std::cos(theta) * std::sin(phi) / r * spherical[1]
      + std::cos(phi) / r / std::sin(theta) * spherical[2];
    cartesian[2] = std::cos(theta) * spherical[0]
      - std::sin(theta) / r * spherical[1];
  }

  void evalMultipole(real_t rho, real_t alpha, real_t beta, complex_t * Ynm, complex_t * YnmTheta) {
    real_t x = std::cos(alpha);
    real_t y = std::sin(alpha);
    real_t invY = y == 0 ? 0 : 1 / y;
    real_t fact = 1;
    real_t pn = 1;
    real_t rhom = 1;
    complex_t ei = std::exp(I * beta);
    complex_t eim = 1.0;
    for (int m=0; m<P; m++) {
      real_t p = pn;
      int npn = m * m + 2 * m;
      int nmn = m * m;
      Ynm[npn] = rhom * p * eim;
      Ynm[nmn] = std::conj(Ynm[npn]);
      real_t p1 = p;
      p = x * (2 * m + 1) * p1;
      YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) * invY * eim;
      rhom *= rho;
      real_t rhon = rhom;
      for (int n=m+1; n<P; n++) {
        int npm = n * n + n + m;
        int nmm = n * n + n - m;
        rhon /= -(n + m);
        Ynm[npm] = rhon * p * eim;
        Ynm[nmm] = std::conj(Ynm[npm]);
        real_t p2 = p1;
        p1 = p;
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
        YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) * invY * eim;
        rhon *= rho;
      }
      rhom /= -(2 * m + 2) * (2 * m + 1);
      pn = -pn * fact * y;
      fact += 2;
      eim *= ei;
    }
  }

  void evalLocal(real_t rho, real_t alpha, real_t beta, complex_t * Ynm) {
    real_t x = std::cos(alpha);
    real_t y = std::sin(alpha);
    real_t fact = 1;
    real_t pn = 1;
    real_t invR = -1.0 / rho;
    real_t rhom = -invR;
    complex_t ei = std::exp(I * beta);
    complex_t eim = 1.0;
    for (int m=0; m<P; m++) {
      real_t p = pn;
      int npn = m * m + 2 * m;
      int nmn = m * m;
      Ynm[npn] = rhom * p * eim;
      Ynm[nmn] = std::conj(Ynm[npn]);
      real_t p1 = p;
      p = x * (2 * m + 1) * p1;
      rhom *= invR;
      real_t rhon = rhom;
      for (int n=m+1; n<P; n++) {
        int npm = n * n + n + m;
        int nmm = n * n + n - m;
        Ynm[npm] = rhon * p * eim;
        Ynm[nmm] = std::conj(Ynm[npm]);
        real_t p2 = p1;
        p1 = p;
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
        rhon *= invR * (n - m + 1);
      }
      pn = -pn * fact * y;
      fact += 2;
      eim *= ei;
    }
  }

  void initKernel() {
    NTERM = P * (P + 1) / 2;
  }

  void P2P(Cell * Ci, Cell * Cj) {
    Body * Bi = Ci->body;
    Body * Bj = Cj->body;
    for (int i=0; i<Ci->numBodies; i++) {
      real_t p = 0;
      vec3 F = 0;
      for (int j=0; j<Cj->numBodies; j++) {
        vec3 dX;
        for (int d=0; d<3; d++) dX[d] = Bi[i].X[d] - Bj[j].X[d] - IX[d] * CYCLE;
        real_t R2 = norm(dX);
        if (R2 != 0) {
          real_t invR2 = 1.0 / R2;
          real_t invR = Bj[j].q * sqrt(invR2);
          p += invR;
          F += dX * invR2 * invR;
        }
      }
#pragma omp atomic
      Bi[i].p += p;
      for (int d=0; d<3; d++) {
#pragma omp atomic
        Bi[i].F[d] -= F[d];
      }
    }
  }

  void P2M(Cell * C) {
    complex_t Ynm[P*P], YnmTheta[P*P];
    for (Body * B=C->body; B!=C->body+C->numBodies; B++) {
      vec3 dX = B->X - C->X;
      real_t rho, alpha, beta;
      cart2sph(dX, rho, alpha, beta);
      evalMultipole(rho, alpha, -beta, Ynm, YnmTheta);
      for (int n=0; n<P; n++) {
        for (int m=0; m<=n; m++) {
          int nm  = n * n + n + m;
          int nms = n * (n + 1) / 2 + m;
          C->M[nms] += B->q * Ynm[nm];
        }
      }
    }
  }

  void M2M(Cell * Ci) {
    complex_t Ynm[P*P], YnmTheta[P*P];
    for (Cell * Cj=Ci->child; Cj!=Ci->child+Ci->numChilds; Cj++) {
      vec3 dX = Ci->X - Cj->X;
      real_t rho, alpha, beta;
      cart2sph(dX, rho, alpha, beta);
      evalMultipole(rho, alpha, beta, Ynm, YnmTheta);
      for (int j=0; j<P; j++) {
        for (int k=0; k<=j; k++) {
          int jks = j * (j + 1) / 2 + k;
          complex_t M = 0;
          for (int n=0; n<=j; n++) {
            for (int m=std::max(-n,-j+k+n); m<=std::min(k-1,n); m++) {
              int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
              int nm    = n * n + n - m;
              M += Cj->M[jnkms] * Ynm[nm] * real_t(ipow2n(m) * oddOrEven(n));
            }
            for (int m=k; m<=std::min(n,j+k-n); m++) {
              int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
              int nm    = n * n + n - m;
              M += std::conj(Cj->M[jnkms]) * Ynm[nm] * real_t(oddOrEven(k+n+m));
            }
          }
          Ci->M[jks] += M;
        }
      }
    }
  }

  void M2L(Cell * Ci, Cell * Cj) {
    complex_t Ynm2[4*P*P];
    vec3 dX;
    for (int d=0; d<3; d++) dX[d] = Ci->X[d] - Cj->X[d] - IX[d] * CYCLE;
    real_t rho, alpha, beta;
    cart2sph(dX, rho, alpha, beta);
    evalLocal(rho, alpha, beta, Ynm2);
    for (int j=0; j<P; j++) {
      real_t Cnm = oddOrEven(j);
      for (int k=0; k<=j; k++) {
        int jks = j * (j + 1) / 2 + k;
        complex_t L = 0;
        for (int n=0; n<P; n++) {
          for (int m=-n; m<0; m++) {
            int nms  = n * (n + 1) / 2 - m;
            int jnkm = (j + n) * (j + n) + j + n + m - k;
            L += std::conj(Cj->M[nms]) * Cnm * Ynm2[jnkm];
          }
          for (int m=0; m<=n; m++) {
            int nms  = n * (n + 1) / 2 + m;
            int jnkm = (j + n) * (j + n) + j + n + m - k;
            real_t Cnm2 = Cnm * oddOrEven((k-m)*(k<m)+m);
            L += Cj->M[nms] * Cnm2 * Ynm2[jnkm];
          }
        }
        Ci->L[jks] += L;
      }
    }
  }

  void L2L(Cell * Cj) {
    complex_t Ynm[P*P], YnmTheta[P*P];
    for (Cell * Ci=Cj->child; Ci!=Cj->child+Cj->numChilds; Ci++) {
      vec3 dX = Ci->X - Cj->X;
      real_t rho, alpha, beta;
      cart2sph(dX, rho, alpha, beta);
      evalMultipole(rho, alpha, beta, Ynm, YnmTheta);
      for (int j=0; j<P; j++) {
        for (int k=0; k<=j; k++) {
          int jks = j * (j + 1) / 2 + k;
          complex_t L = 0;
          for (int n=j; n<P; n++) {
            for (int m=j+k-n; m<0; m++) {
              int jnkm = (n - j) * (n - j) + n - j + m - k;
              int nms  = n * (n + 1) / 2 - m;
              L += std::conj(Cj->L[nms]) * Ynm[jnkm] * real_t(oddOrEven(k));
            }
            for (int m=0; m<=n; m++) {
              if (n-j >= abs(m-k)) {
                int jnkm = (n - j) * (n - j) + n - j + m - k;
                int nms  = n * (n + 1) / 2 + m;
                L += Cj->L[nms] * Ynm[jnkm] * real_t(oddOrEven((m-k)*(m<k)));
              }
            }
          }
          Ci->L[jks] += L;
        }
      }
    }
  }

  void L2P(Cell * Ci) {
    complex_t Ynm[P*P], YnmTheta[P*P];
    for (Body * B=Ci->body; B!=Ci->body+Ci->numBodies; B++) {
      vec3 dX = B->X - Ci->X;
      vec3 spherical = 0;
      vec3 cartesian = 0;
      real_t r, theta, phi;
      cart2sph(dX, r, theta, phi);
      evalMultipole(r, theta, phi, Ynm, YnmTheta);
      for (int n=0; n<P; n++) {
        int nm  = n * n + n;
        int nms = n * (n + 1) / 2;
        B->p += std::real(Ci->L[nms] * Ynm[nm]);
        spherical[0] += std::real(Ci->L[nms] * Ynm[nm]) / r * n;
        spherical[1] += std::real(Ci->L[nms] * YnmTheta[nm]);
        for (int m=1; m<=n; m++) {
          nm  = n * n + n + m;
          nms = n * (n + 1) / 2 + m;
          B->p += 2 * std::real(Ci->L[nms] * Ynm[nm]);
          spherical[0] += 2 * std::real(Ci->L[nms] * Ynm[nm]) / r * n;
          spherical[1] += 2 * std::real(Ci->L[nms] * YnmTheta[nm]);
          spherical[2] += 2 * std::real(Ci->L[nms] * Ynm[nm] * I) * m;
        }
      }
      sph2cart(r, theta, phi, spherical, cartesian);
      B->F += cartesian;
    }
  }
}
#endif
