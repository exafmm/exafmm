#ifndef laplace_h
#define laplace_h
#include "namespace.h"
#include "types.h"

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
    void cart2sph(vec3 dX, real_t & r, real_t & theta, real_t & phi);
    //! Spherical to cartesian coordinates
    void sph2cart(real_t r, real_t theta, real_t phi, vec3 spherical, vec3 & cartesian);
    //! Evaluate solid harmonics \f$ r^n Y_{n}^{m} \f$
    void evalMultipole(real_t rho, real_t alpha, real_t beta, complex_t * Ynm, complex_t * YnmTheta);
    //! Evaluate singular harmonics \f$ r^{-n-1} Y_n^m \f$
    void evalLocal(real_t rho, real_t alpha, real_t beta, complex_t * Ynm2);

  public:
    Kernel(int _P, real_t _eps2, complex_t _wavek);
    void P2P(Target* Bi, int ni, Source* Bj, int nj);
    void P2M(vec3 & X, Coefs & M, Source* B, int nj);
    void M2M(vec3 & Xi, Coefs & Mi, vec3 & Xj, Coefs & Mj);
    void M2L(vec3 & Xi, Coefs & Li, vec3 & Xj, Coefs & Mj);
    void L2L(vec3 & Xi, Coefs & Li, vec3 & Xj, Coefs & Lj);
    void L2P(Target* B, int ni, vec3 & X, Coefs & L);
  };
}
#endif
