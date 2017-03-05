#ifndef laplace_h
#define laplace_h
#include "namespace.h"

namespace EXAFMM_NAMESPACE {
  class Kernel {
  private:
    std::vector<real_t> prefactor;                              // sqrt( (n - |m|)! / (n + |m|)! )
    std::vector<real_t> Anm;                                    // (-1)^n / sqrt( (n + m)! / (n - m)! )
    std::vector<complex_t> Cnm;                                 // M2L translation matrix Cjknm

  public:
    const int P;
    const int NTERM;
    real_t eps2;
    complex_t wavek;
    vec3 Xperiodic;

  private:
    //! Odd or even
    inline int oddOrEven(int n) {
      return (((n) & 1) == 1) ? -1 : 1;
    }

    //! Get r,theta,phi from x,y,z
    void cart2sph(vec3 dX, real_t & r, real_t & theta, real_t & phi) {
      r = sqrt(norm(dX));                                       // r = sqrt(x^2 + y^2 + z^2)
      theta = r == 0 ? 0 : acos(dX[2] / r);                     // theta = acos(z / r)
      phi = atan2(dX[1], dX[0]);                                // phi = atan(y / x)
    }

    //! Spherical to cartesian coordinates
    void sph2cart(real_t r, real_t theta, real_t phi, vec3 spherical, vec3 & cartesian) {
      cartesian[0] = std::sin(theta) * std::cos(phi) * spherical[0] // x component (not x itself)
        + std::cos(theta) * std::cos(phi) / r * spherical[1]
        - std::sin(phi) / r / std::sin(theta) * spherical[2];
      cartesian[1] = std::sin(theta) * std::sin(phi) * spherical[0] // y component (not y itself)
        + std::cos(theta) * std::sin(phi) / r * spherical[1]
        + std::cos(phi) / r / std::sin(theta) * spherical[2];
      cartesian[2] = std::cos(theta) * spherical[0]             // z component (not z itself)
        - std::sin(theta) / r * spherical[1];
    }

    //! Evaluate solid harmonics \f$ r^n Y_{n}^{m} \f$
    void evalMultipole(real_t rho, real_t alpha, real_t beta, complex_t * Ynm, complex_t * YnmTheta) {
      real_t x = std::cos(alpha);                               // x = cos(alpha)
      real_t y = std::sin(alpha);                               // y = sin(alpha)
      real_t fact = 1;                                          // Initialize 2 * m + 1
      real_t pn = 1;                                            // Initialize Legendre polynomial Pn
      real_t rhom = 1;                                          // Initialize rho^m
      for (int m=0; m<P; m++) {                                 // Loop over m in Ynm
        complex_t eim = std::exp(I * real_t(m * beta));         //  exp(i * m * beta)
        real_t p = pn;                                          //  Associated Legendre polynomial Pnm
        int npn = m * m + 2 * m;                                //  Index of Ynm for m > 0
        int nmn = m * m;                                        //  Index of Ynm for m < 0
        Ynm[npn] = rhom * p * prefactor[npn] * eim;             //  rho^m * Ynm for m > 0
        Ynm[nmn] = std::conj(Ynm[npn]);                         //  Use conjugate relation for m < 0
        real_t p1 = p;                                          //  Pnm-1
        p = x * (2 * m + 1) * p1;                               //  Pnm using recurrence relation
        YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) / y * prefactor[npn] * eim;// theta derivative of r^n * Ynm
        rhom *= rho;                                            //  rho^m
        real_t rhon = rhom;                                     //  rho^n
        for (int n=m+1; n<P; n++) {                             //  Loop over n in Ynm
          int npm = n * n + n + m;                              //   Index of Ynm for m > 0
          int nmm = n * n + n - m;                              //   Index of Ynm for m < 0
          Ynm[npm] = rhon * p * prefactor[npm] * eim;           //   rho^n * Ynm
          Ynm[nmm] = std::conj(Ynm[npm]);                       //   Use conjugate relation for m < 0
          real_t p2 = p1;                                       //   Pnm-2
          p1 = p;                                               //   Pnm-1
          p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);//   Pnm using recurrence relation
          YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) / y * prefactor[npm] * eim;// theta derivative
          rhon *= rho;                                          //   Update rho^n
        }                                                       //  End loop over n in Ynm
        pn = -pn * fact * y;                                    //  Pn
        fact += 2;                                              //  2 * m + 1
      }                                                         // End loop over m in Ynm
    }

  //! Evaluate singular harmonics \f$ r^{-n-1} Y_n^m \f$
    void evalLocal(real_t rho, real_t alpha, real_t beta, complex_t * Ynm2) {
      real_t x = std::cos(alpha);                               // x = cos(alpha)
      real_t y = std::sin(alpha);                               // y = sin(alpha)
      real_t fact = 1;                                          // Initialize 2 * m + 1
      real_t pn = 1;                                            // Initialize Legendre polynomial Pn
      real_t rhom = 1.0 / rho;                                  // Initialize rho^(-m-1)
      for (int m=0; m<2*P; m++) {                               // Loop over m in Ynm
        complex_t eim = std::exp(I * real_t(m * beta));         //  exp(i * m * beta)
        real_t p = pn;                                          //  Associated Legendre polynomial Pnm
        int npn = m * m + 2 * m;                                //  Index of Ynm for m > 0
        int nmn = m * m;                                        //  Index of Ynm for m < 0
        Ynm2[npn] = rhom * p * prefactor[npn] * eim;            //  rho^(-m-1) * Ynm for m > 0
        Ynm2[nmn] = std::conj(Ynm2[npn]);                       //  Use conjugate relation for m < 0
        real_t p1 = p;                                          //  Pnm-1
        p = x * (2 * m + 1) * p1;                               //  Pnm using recurrence relation
        rhom /= rho;                                            //  rho^(-m-1)
        real_t rhon = rhom;                                     //  rho^(-n-1)
        for (int n=m+1; n<2*P; n++) {                           //  Loop over n in Ynm
          int npm = n * n + n + m;                              //   Index of Ynm for m > 0
          int nmm = n * n + n - m;                              //   Index of Ynm for m < 0
          Ynm2[npm] = rhon * p * prefactor[npm] * eim;          //   rho^n * Ynm for m > 0
          Ynm2[nmm] = std::conj(Ynm2[npm]);                     //   Use conjugate relation for m < 0
          real_t p2 = p1;                                       //   Pnm-2
          p1 = p;                                               //   Pnm-1
          p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);//   Pnm using recurrence relation
          rhon /= rho;                                          //   rho^(-n-1)
        }                                                       //  End loop over n in Ynm
        pn = -pn * fact * y;                                    //  Pn
        fact += 2;                                              //  2 * m + 1
      }                                                         // End loop over m in Ynm
    }

  public:
    Kernel(int _P, real_t _eps2, complex_t _wavek) : P(_P), NTERM(P*(P+1)/2), eps2(_eps2), wavek(_wavek) {
      Xperiodic = 0;
      prefactor.resize(4*P*P);
      Anm.resize(4*P*P);
      Cnm.resize(P*P*P*P);
      for (int n=0; n<2*P; n++) {                               // Loop over n in Anm
        for (int m=-n; m<=n; m++) {                             //  Loop over m in Anm
          int nm = n*n+n+m;                                     //   Index of Anm
          int nabsm = abs(m);                                   //   |m|
          real_t fnmm = EPS;                                    //   Initialize (n - m)!
          for (int i=1; i<=n-m; i++) fnmm *= i;                 //   (n - m)!
          real_t fnpm = EPS;                                    //   Initialize (n + m)!
          for (int i=1; i<=n+m; i++) fnpm *= i;                 //   (n + m)!
          real_t fnma = 1.0;                                    //   Initialize (n - |m|)!
          for (int i=1; i<=n-nabsm; i++) fnma *= i;             //   (n - |m|)!
          real_t fnpa = 1.0;                                    //   Initialize (n + |m|)!
          for (int i=1; i<=n+nabsm; i++) fnpa *= i;             //   (n + |m|)!
          prefactor[nm] = std::sqrt(fnma/fnpa);                 //   sqrt( (n - |m|)! / (n + |m|)! )
          Anm[nm] = oddOrEven(n)/std::sqrt(fnmm*fnpm);          //   (-1)^n / sqrt( (n + m)! / (n - m)! )
        }                                                       //  End loop over m in Anm
      }                                                         // End loop over n in Anm
      for (int j=0, jk=0, jknm=0; j<P; j++) {                   // Loop over j in Cjknm
        for (int k=-j; k<=j; k++, jk++) {                       //  Loop over k in Cjknm
          for (int n=0, nm=0; n<P; n++) {                       //   Loop over n in Cjknm
            for (int m=-n; m<=n; m++, nm++, jknm++) {           //    Loop over m in Cjknm
              const int jnkm = (j+n)*(j+n)+j+n+m-k;             //     Index C_{j+n}^{m-k}
              Cnm[jknm] = std::pow(I,real_t(abs(k-m)-abs(k)-abs(m)))//     Cjknm
                * real_t(oddOrEven(j)*Anm[nm]*Anm[jk]/Anm[jnkm]) * EPS;
            }                                                   //    End loop over m in Cjknm
          }                                                     //   End loop over n in Cjknm
        }                                                       //  End loop over in k in Cjknm
      }                                                         // End loop over in j in Cjknm
    }

    void P2P(C_iter Ci, C_iter Cj) {
      B_iter Bi = Ci->BODY;
      B_iter Bj = Cj->BODY;
      int ni = Ci->NBODY;
      int nj = Cj->NBODY;
      for (int i=0; i<ni; i++) {
        kreal_t pot = 0;
        kreal_t ax = 0;
        kreal_t ay = 0;
        kreal_t az = 0;
        for (int j=0; j<nj; j++) {
          vec3 dX = Bi[i].X - Bj[j].X - Xperiodic;
          real_t R2 = norm(dX) + eps2;
          if (R2 != 0) {
            real_t invR2 = 1.0 / R2;
            real_t invR = Bj[j].SRC * sqrt(invR2);
            dX *= invR2 * invR;
            pot += invR;
            ax += dX[0];
            ay += dX[1];
            az += dX[2];
          }
        }
        Bi[i].TRG[0] += pot;
        Bi[i].TRG[1] -= ax;
        Bi[i].TRG[2] -= ay;
        Bi[i].TRG[3] -= az;
      }
    }

    void P2M(C_iter C) {
      complex_t Ynm[P*P], YnmTheta[P*P];
      for (B_iter B=C->BODY; B!=C->BODY+C->NBODY; B++) {
        vec3 dX = B->X - C->X;
        real_t rho, alpha, beta;
        cart2sph(dX, rho, alpha, beta);
        evalMultipole(rho, alpha, -beta, Ynm, YnmTheta);
        for (int n=0; n<P; n++) {
          for (int m=0; m<=n; m++) {
            int nm  = n * n + n + m;
            int nms = n * (n + 1) / 2 + m;
            C->M[nms] += B->SRC * Ynm[nm];
          }
        }
      }
    }

    void M2M(C_iter Ci, C_iter C0) {
      complex_t Ynm[P*P], YnmTheta[P*P];
      for (C_iter Cj=C0+Ci->ICHILD; Cj!=C0+Ci->ICHILD+Ci->NCHILD; Cj++) {
        vec3 dX = Ci->X - Cj->X;
        real_t rho, alpha, beta;
        cart2sph(dX, rho, alpha, beta);
        evalMultipole(rho, alpha, -beta, Ynm, YnmTheta);
        for (int j=0; j<P; j++) {
          for (int k=0; k<=j; k++) {
            int jk = j * j + j + k;
            int jks = j * (j + 1) / 2 + k;
            complex_t M = 0;
            for (int n=0; n<=j; n++) {
              for (int m=-n; m<=std::min(k-1,n); m++) {
                if (j-n >= k-m) {
                  int jnkm  = (j - n) * (j - n) + j - n + k - m;
                  int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
                  int nm    = n * n + n + m;
                  M += Cj->M[jnkms] * std::pow(I,real_t(m-abs(m))) * Ynm[nm]
                    * real_t(oddOrEven(n) * Anm[nm] * Anm[jnkm] / Anm[jk]);
                }
              }
              for (int m=k; m<=n; m++) {
                if (j-n >= m-k) {
                  int jnkm  = (j - n) * (j - n) + j - n + k - m;
                  int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
                  int nm    = n * n + n + m;
                  M += std::conj(Cj->M[jnkms]) * Ynm[nm]
                    * real_t(oddOrEven(k+n+m) * Anm[nm] * Anm[jnkm] / Anm[jk]);
                }
              }
            }
            Ci->M[jks] += M * EPS;
          }
        }
      }
    }

    void M2L(C_iter Ci, C_iter Cj) {
      complex_t Ynm2[4*P*P];
      vec3 dX = Ci->X - Cj->X - Xperiodic;
      real_t rho, alpha, beta;
      cart2sph(dX, rho, alpha, beta);
      evalLocal(rho, alpha, beta, Ynm2);
      for (int j=0; j<P; j++) {
        for (int k=0; k<=j; k++) {
          int jk = j * j + j + k;
          int jks = j * (j + 1) / 2 + k;
          complex_t L = 0;
          for (int n=0; n<P; n++) {
            for (int m=-n; m<0; m++) {
              int nm   = n * n + n + m;
              int nms  = n * (n + 1) / 2 - m;
              int jknm = jk * P * P + nm;
              int jnkm = (j + n) * (j + n) + j + n + m - k;
              L += std::conj(Cj->M[nms]) * Cnm[jknm] * Ynm2[jnkm];
            }
            for (int m=0; m<=n; m++) {
              int nm   = n * n + n + m;
              int nms  = n * (n + 1) / 2 + m;
              int jknm = jk * P * P + nm;
              int jnkm = (j + n) * (j + n) + j + n + m - k;
              L += Cj->M[nms] * Cnm[jknm] * Ynm2[jnkm];
            }
          }
          Ci->L[jks] += L;
        }
      }
    }

    void L2L(C_iter Ci, C_iter C0) {
      complex_t Ynm[P*P], YnmTheta[P*P];
      C_iter Cj = C0 + Ci->IPARENT;
      vec3 dX = Ci->X - Cj->X;
      real_t rho, alpha, beta;
      cart2sph(dX, rho, alpha, beta);
      evalMultipole(rho, alpha, beta, Ynm, YnmTheta);
      for (int j=0; j<P; j++) {
        for (int k=0; k<=j; k++) {
          int jk = j * j + j + k;
          int jks = j * (j + 1) / 2 + k;
          complex_t L = 0;
          for (int n=j; n<P; n++) {
            for (int m=j+k-n; m<0; m++) {
              int jnkm = (n - j) * (n - j) + n - j + m - k;
              int nm   = n * n + n - m;
              int nms  = n * (n + 1) / 2 - m;
              L += std::conj(Cj->L[nms]) * Ynm[jnkm]
                * real_t(oddOrEven(k) * Anm[jnkm] * Anm[jk] / Anm[nm]);
            }
            for (int m=0; m<=n; m++) {
              if (n-j >= abs(m-k)) {
                int jnkm = (n - j) * (n - j) + n - j + m - k;
                int nm   = n * n + n + m;
                int nms  = n * (n + 1) / 2 + m;
                L += Cj->L[nms] * std::pow(I,real_t(m-k-abs(m-k)))
                  * Ynm[jnkm] * Anm[jnkm] * Anm[jk] / Anm[nm];
              }
            }
          }
          Ci->L[jks] += L * EPS;
        }
      }
    }

    void L2P(C_iter Ci) {
      complex_t Ynm[P*P], YnmTheta[P*P];
      for (B_iter B=Ci->BODY; B!=Ci->BODY+Ci->NBODY; B++) {
        vec3 dX = B->X - Ci->X + EPS;
        vec3 spherical = 0;
        vec3 cartesian = 0;
        real_t r, theta, phi;
        cart2sph(dX, r, theta, phi);
        evalMultipole(r, theta, phi, Ynm, YnmTheta);
        for (int n=0; n<P; n++) {
          int nm  = n * n + n;
          int nms = n * (n + 1) / 2;
          B->TRG[0] += std::real(Ci->L[nms] * Ynm[nm]);
          spherical[0] += std::real(Ci->L[nms] * Ynm[nm]) / r * n;
          spherical[1] += std::real(Ci->L[nms] * YnmTheta[nm]);
          for (int m=1; m<=n; m++) {
            nm  = n * n + n + m;
            nms = n * (n + 1) / 2 + m;
            B->TRG[0] += 2 * std::real(Ci->L[nms] * Ynm[nm]);
            spherical[0] += 2 * std::real(Ci->L[nms] * Ynm[nm]) / r * n;
            spherical[1] += 2 * std::real(Ci->L[nms] * YnmTheta[nm]);
            spherical[2] += 2 * std::real(Ci->L[nms] * Ynm[nm] * I) * m;
          }
        }
        sph2cart(r, theta, phi, spherical, cartesian);
        B->TRG[1] += cartesian[0];
        B->TRG[2] += cartesian[1];
        B->TRG[3] += cartesian[2];
      }
    }
  };
}
#endif
