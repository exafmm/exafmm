#ifndef kernel_h
#define kernel_h
#include "exafmm.h"

namespace exafmm {
  const complex_t I(0.,1.);                                     //!< Imaginary unit

  //!< L2 norm of vector X
  inline real_t norm(real_t * X) {
    return X[0] * X[0] + X[1] * X[1] + X[2] * X[2];             // L2 norm
  }

  //! Odd or even
  inline int oddOrEven(int n) {
    return (((n) & 1) == 1) ? -1 : 1;                           // Odd: -1, Even: 1
  }

  //! i^2n
  inline int ipow2n(int n) {
    return (n >= 0) ? 1 : oddOrEven(n);                         // i^2n
  }

  //! Get r,theta,phi from x,y,z
  void cart2sph(real_t * dX, real_t & r, real_t & theta, real_t & phi) {
    r = sqrt(norm(dX));                                         // r = sqrt(x^2 + y^2 + z^2)
    theta = r == 0 ? 0 : acos(dX[2] / r);                       // theta = acos(z / r)
    phi = atan2(dX[1], dX[0]);                                  // phi = atan(y / x)
  }

  //! Spherical to cartesian coordinates
  void sph2cart(real_t r, real_t theta, real_t phi, real_t * spherical, real_t * cartesian) {
    cartesian[0] = std::sin(theta) * std::cos(phi) * spherical[0]// x component (not x itself)
      + std::cos(theta) * std::cos(phi) / r * spherical[1]
      - std::sin(phi) / r / std::sin(theta) * spherical[2];
    cartesian[1] = std::sin(theta) * std::sin(phi) * spherical[0]// y component (not y itself)
      + std::cos(theta) * std::sin(phi) / r * spherical[1]
      + std::cos(phi) / r / std::sin(theta) * spherical[2];
    cartesian[2] = std::cos(theta) * spherical[0]               // z component (not z itself)
      - std::sin(theta) / r * spherical[1];
  }

  //! Evaluate solid harmonics \f$ r^n Y_{n}^{m} \f$
  void evalMultipole(real_t rho, real_t alpha, real_t beta, complex_t * Ynm, complex_t * YnmTheta) {
    real_t x = std::cos(alpha);                                 // x = cos(alpha)
    real_t y = std::sin(alpha);                                 // y = sin(alpha)
    real_t invY = y == 0 ? 0 : 1 / y;                           // 1 / y
    real_t fact = 1;                                            // Initialize 2 * m + 1
    real_t pn = 1;                                              // Initialize Legendre polynomial Pn
    real_t rhom = 1;                                            // Initialize rho^m
    complex_t ei = std::exp(I * beta);                          // exp(i * beta)
    complex_t eim = 1.0;                                        // Initialize exp(i * m * beta)
    for (int m=0; m<P; m++) {                                   // Loop over m in Ynm
      real_t p = pn;                                            //  Associated Legendre polynomial Pnm
      int npn = m * m + 2 * m;                                  //  Index of Ynm for m > 0
      int nmn = m * m;                                          //  Index of Ynm for m < 0
      Ynm[npn] = rhom * p * eim;                                //  rho^m * Ynm for m > 0
      Ynm[nmn] = std::conj(Ynm[npn]);                           //  Use conjugate relation for m < 0
      real_t p1 = p;                                            //  Pnm-1
      p = x * (2 * m + 1) * p1;                                 //  Pnm using recurrence relation
      YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) * invY * eim; //  theta derivative of r^n * Ynm
      rhom *= rho;                                              //  rho^m
      real_t rhon = rhom;                                       //  rho^n
      for (int n=m+1; n<P; n++) {                               //  Loop over n in Ynm
        int npm = n * n + n + m;                                //   Index of Ynm for m > 0
        int nmm = n * n + n - m;                                //   Index of Ynm for m < 0
        rhon /= -(n + m);                                       //   Update factorial
        Ynm[npm] = rhon * p * eim;                              //   rho^n * Ynm
        Ynm[nmm] = std::conj(Ynm[npm]);                         //   Use conjugate relation for m < 0
        real_t p2 = p1;                                         //   Pnm-2
        p1 = p;                                                 //   Pnm-1
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);//   Pnm using recurrence relation
        YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) * invY * eim;// theta derivative
        rhon *= rho;                                            //   Update rho^n
      }                                                         //  End loop over n in Ynm
      rhom /= -(2 * m + 2) * (2 * m + 1);                       //  Update factorial
      pn = -pn * fact * y;                                      //  Pn
      fact += 2;                                                //  2 * m + 1
      eim *= ei;                                                //  Update exp(i * m * beta)
    }                                                           // End loop over m in Ynm
  }

  //! Evaluate singular harmonics \f$ r^{-n-1} Y_n^m \f$
  void evalLocal(real_t rho, real_t alpha, real_t beta, complex_t * Ynm) {
    real_t x = std::cos(alpha);                                 // x = cos(alpha)
    real_t y = std::sin(alpha);                                 // y = sin(alpha)
    real_t fact = 1;                                            // Initialize 2 * m + 1
    real_t pn = 1;                                              // Initialize Legendre polynomial Pn
    real_t invR = -1.0 / rho;                                   // - 1 / rho
    real_t rhom = -invR;                                        // Initialize rho^(-m-1)
    complex_t ei = std::exp(I * beta);                          // exp(i * beta)
    complex_t eim = 1.0;                                        // Initialize exp(i * m * beta)
    for (int m=0; m<P; m++) {                                   // Loop over m in Ynm
      real_t p = pn;                                            //  Associated Legendre polynomial Pnm
      int npn = m * m + 2 * m;                                  //  Index of Ynm for m > 0
      int nmn = m * m;                                          //  Index of Ynm for m < 0
      Ynm[npn] = rhom * p * eim;                                //  rho^(-m-1) * Ynm for m > 0
      Ynm[nmn] = std::conj(Ynm[npn]);                           //  Use conjugate relation for m < 0
      real_t p1 = p;                                            //  Pnm-1
      p = x * (2 * m + 1) * p1;                                 //  Pnm using recurrence relation
      rhom *= invR;                                             //  rho^(-m-1)
      real_t rhon = rhom;                                       //  rho^(-n-1)
      for (int n=m+1; n<P; n++) {                               //  Loop over n in Ynm
        int npm = n * n + n + m;                                //   Index of Ynm for m > 0
        int nmm = n * n + n - m;                                //   Index of Ynm for m < 0
        Ynm[npm] = rhon * p * eim;                              //   rho^n * Ynm for m > 0
        Ynm[nmm] = std::conj(Ynm[npm]);                         //   Use conjugate relation for m < 0
        real_t p2 = p1;                                         //   Pnm-2
        p1 = p;                                                 //   Pnm-1
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);//   Pnm using recurrence relation
        rhon *= invR * (n - m + 1);                             //   rho^(-n-1)
      }                                                         //  End loop over n in Ynm
      pn = -pn * fact * y;                                      //  Pn
      fact += 2;                                                //  2 * m + 1
      eim *= ei;                                                //  Update exp(i * m * beta)
    }                                                           // End loop over m in Ynm
  }

  void initKernel() {
    NTERM = P * (P + 1) / 2;                                    // Calculate number of coefficients
  }

  void P2P(Cell * Ci, Cell * Cj) {
    Body * Bi = Ci->BODY;
    Body * Bj = Cj->BODY;
    int ni = Ci->NBODY;
    int nj = Cj->NBODY;
    for (int i=0; i<ni; i++) {
      real_t pot = 0;
      real_t ax = 0;
      real_t ay = 0;
      real_t az = 0;
      for (int j=0; j<nj; j++) {
        for (int d=0; d<3; d++) dX[d] = Bi[i].X[d] - Bj[j].X[d] - iX[d] * cycle;
        real_t R2 = norm(dX);
        if (R2 != 0) {
          real_t invR2 = 1.0 / R2;
          real_t invR = Bj[j].q * sqrt(invR2);
          for (int d=0; d<3; d++) dX[d] *= invR2 * invR;
          pot += invR;
          ax += dX[0];
          ay += dX[1];
          az += dX[2];
        }
      }
      Bi[i].p += pot;
      Bi[i].F[0] -= ax;
      Bi[i].F[1] -= ay;
      Bi[i].F[2] -= az;
    }
  }

  void P2M(Cell * C) {
    complex_t Ynm[P*P], YnmTheta[P*P];
    for (Body * B=C->BODY; B!=C->BODY+C->NBODY; B++) {
      for (int d=0; d<3; d++) dX[d] = B->X[d] - C->X[d];
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
    for (Cell * Cj=Ci->CHILD; Cj!=Ci->CHILD+Ci->NCHILD; Cj++) {
      for (int d=0; d<3; d++) dX[d] = Ci->X[d] - Cj->X[d];
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
    for (int d=0; d<3; d++) dX[d] = Ci->X[d] - Cj->X[d] - iX[d] * cycle;
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
    for (Cell * Ci=Cj->CHILD; Ci!=Cj->CHILD+Cj->NCHILD; Ci++) {
      for (int d=0; d<3; d++) dX[d] = Ci->X[d] - Cj->X[d];
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
    for (Body * B=Ci->BODY; B!=Ci->BODY+Ci->NBODY; B++) {
      for (int d=0; d<3; d++) dX[d] = B->X[d] - Ci->X[d];
      real_t spherical[3] = {0, 0, 0};
      real_t cartesian[3] = {0, 0, 0};
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
      B->F[0] += cartesian[0];
      B->F[1] += cartesian[1];
      B->F[2] += cartesian[2];
    }
  }
}
#endif
