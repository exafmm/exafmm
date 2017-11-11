#ifndef stokes_h
#define stokes_h
#include "exafmm.h"

namespace exafmm {
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
    NTERM = 4 * P * (P + 1) / 2;
  }

  void P2P(Cell * Ci, Cell * Cj) {
    Body * Bi = Ci->body;
    Body * Bj = Cj->body;
    for (int i=0; i<Ci->numBodies; i++) {
      vec3 F = 0;
      for (int j=0; j<Cj->numBodies; j++) {
        vec3 dX = Bi[i].X - Bj[j].X;
        real_t R2 = norm(dX);
        if (R2 != 0) {
          real_t fdX = Bj[j].q[0] * dX[0] + Bj[j].q[1] * dX[1] + Bj[j].q[2] * dX[2];
          real_t invR = 1.0 / std::sqrt(R2);
          real_t invR3 = invR / R2;
          F += Bj[j].q * invR;
          F += dX * invR3 * fdX;
        }
      }
      Bi[i].F += F;
    }
  }
  
  void P2M(Cell * C) {
    complex_t Ynm[P*P], YnmTheta[P*P];
    for (Body * B=C->body; B!=C->body+C->numBodies; B++) {
      vec3 dX = B->X - C->X;
      real_t rho, alpha, beta;
      cart2sph(dX, rho, alpha, beta);
      evalMultipole(rho, alpha, -beta, Ynm, YnmTheta);
      real_t fX = B->q[0] * B->X[0] + B->q[1] * B->X[1] + B->q[2] * B->X[2];
      for (int n=0; n<P; n++) {
        for (int m=0; m<=n; m++) {
          int nm  = n * n + n + m;
          int nms = n * (n + 1) / 2 + m;
          for (int d=0; d<3; d++) {
            C->M[4*nms+d] += B->q[d] * Ynm[nm];
          }
          C->M[4*nms+3] += fX * Ynm[nm];
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
          vec<4,complex_t> M = complex_t(0.,0.);
          for (int n=0; n<=j; n++) {
            for (int m=std::max(-n,-j+k+n); m<=std::min(k-1,n); m++) {
              int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
              int nm    = n * n + n - m;
              for (int d=0; d<4; d++) {
                M[d] += Cj->M[4*jnkms+d] * Ynm[nm] * real_t(ipow2n(m) * oddOrEven(n));
              }
            }
            for (int m=k; m<=std::min(n,j+k-n); m++) {
              int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
              int nm    = n * n + n - m;
              for (int d=0; d<4; d++) {
                M[d] += std::conj(Cj->M[4*jnkms+d]) * Ynm[nm] * real_t(oddOrEven(k+n+m));
              }
            }
          }
          for (int d=0; d<4; d++) {
            Ci->M[4*jks+d] += M[d];
          }
        }
      }
    }
  }

  void M2L(Cell * Ci, Cell * Cj) {
    complex_t Ynm2[4*P*P];
    vec3 dX = Ci->X - Cj->X;
    real_t rho, alpha, beta;
    cart2sph(dX, rho, alpha, beta);
    evalLocal(rho, alpha, beta, Ynm2);
    for (int j=0; j<P; j++) {
      real_t Cnm = oddOrEven(j);
      for (int k=0; k<=j; k++) {
        int jks = j * (j + 1) / 2 + k;
        vec<4,complex_t> L = complex_t(0.,0.);
        for (int n=0; n<P; n++) {
          for (int m=-n; m<0; m++) {
            int nms  = n * (n + 1) / 2 - m;
            int jnkm = (j + n) * (j + n) + j + n + m - k;
            for (int d=0; d<4; d++) {
              L[d] += std::conj(Cj->M[4*nms+d]) * Cnm * Ynm2[jnkm];
            }
          }
          for (int m=0; m<=n; m++) {
            int nms  = n * (n + 1) / 2 + m;
            int jnkm = (j + n) * (j + n) + j + n + m - k;
            real_t Cnm2 = Cnm * oddOrEven((k-m)*(k<m)+m);
            for (int d=0; d<4; d++) {
              L[d] += Cj->M[4*nms+d] * Cnm2 * Ynm2[jnkm];
            }
          }
        }
        for (int d=0; d<4; d++) {
          Ci->L[4*jks+d] += L[d];
        }
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
          vec<4,complex_t> L = complex_t(0.,0.);
          for (int n=j; n<P; n++) {
            for (int m=j+k-n; m<0; m++) {
              int jnkm = (n - j) * (n - j) + n - j + m - k;
              int nms  = n * (n + 1) / 2 - m;
              for (int d=0; d<4; d++) {
                L[d] += std::conj(Cj->L[4*nms+d]) * Ynm[jnkm] * real_t(oddOrEven(k));
              }
            }
            for (int m=0; m<=n; m++) {
              if (n-j >= abs(m-k)) {
                int jnkm = (n - j) * (n - j) + n - j + m - k;
                int nms  = n * (n + 1) / 2 + m;
                for (int d=0; d<4; d++) {
                  L[d] += Cj->L[4*nms+d] * Ynm[jnkm] * real_t(oddOrEven((m-k)*(m<k)));
                }
              }
            }
          }
          for (int d=0; d<4; d++) {
            Ci->L[4*jks+d] += L[d];
          }
        }
      }
    }
  }

  void L2P(Cell * C) {
    complex_t Ynm[P*P], YnmTheta[P*P];
    for (Body * B=C->body; B!=C->body+C->numBodies; B++) {
      vec3 dX = B->X - C->X;
      vec3 gradient[4] = {0., 0., 0., 0.};
      vec3 cartesian = 0;
      real_t r, theta, phi;
      cart2sph(dX, r, theta, phi);
      evalMultipole(r, theta, phi, Ynm, YnmTheta);
      for (int n=0; n<P; n++) {
        int nm  = n * n + n;
        int nms = n * (n + 1) / 2;
        for (int d=0; d<4; d++) {
          if (d<3) B->F[d] += std::real(C->L[4*nms+d] * Ynm[nm]);
          gradient[d][0] += std::real(C->L[4*nms+d] * Ynm[nm]) / r * n;
          gradient[d][1] += std::real(C->L[4*nms+d] * YnmTheta[nm]);
        }
        for (int m=1; m<=n; m++) {
          nm  = n * n + n + m;
          nms = n * (n + 1) / 2 + m;
          for (int d=0; d<4; d++) {
            if (d<3) B->F[d] += 2 * std::real(C->L[4*nms+d] * Ynm[nm]);
            gradient[d][0] += 2 * std::real(C->L[4*nms+d] * Ynm[nm]) / r * n;
            gradient[d][1] += 2 * std::real(C->L[4*nms+d] * YnmTheta[nm]);
            gradient[d][2] += 2 * std::real(C->L[4*nms+d] * Ynm[nm] * I) * m;
          }
        }
      }
      for (int d=0; d<4; d++) {
        sph2cart(r, theta, phi, gradient[d], cartesian);
        if (d<3) cartesian *= -B->X[d];
        gradient[d] = cartesian;
      }
      B->F[0] += gradient[0][0] + gradient[1][0] + gradient[2][0] + gradient[3][0];
      B->F[1] += gradient[0][1] + gradient[1][1] + gradient[2][1] + gradient[3][1];
      B->F[2] += gradient[0][2] + gradient[1][2] + gradient[2][2] + gradient[3][2];
    }
  }
}
#endif
