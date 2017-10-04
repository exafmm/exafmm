#ifndef helmholtz_h
#define helmholtz_h
#include "exafmm.h"

namespace exafmm {
  int NQUAD, NQUAD2;
  std::vector<real_t> XQUAD, XQUAD2;
  std::vector<real_t> WQUAD, WQUAD2;
  std::vector<real_t> Anm1, Anm2;

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

  void getAnm(std::vector<real_t> & Anm1, std::vector<real_t> & Anm2) {
    Anm1[0] = 1;
    Anm2[0] = 1;
    for (int m=0; m<=P; m++) {
      int ms = m * (m + 1) / 2 + m;
      int mps = (m + 1) * (m + 2) / 2 + m;
      if(m>0) Anm1[ms] = sqrt((2 * m - 1.0) / (2 * m));
      if(m<P) Anm1[mps] = sqrt(2 * m + 1.0);
      for (int n=m+2; n<=P; n++) {
        int nms = n * (n + 1) / 2 + m;
        Anm1[nms] = 2 * n - 1;
        Anm2[nms] = sqrt((n + m - 1.0) * (n - m - 1.0));
        Anm1[nms] /= sqrt(real_t(n - m) * (n + m));
        Anm2[nms] /= sqrt(real_t(n - m) * (n + m));
      }
    }
  }

  void get_Ynm(int nterms, real_t x, real_t * Ynm) {
    real_t y = -sqrt((1 - x) * (1 + x));
    Ynm[0] = 1;
    for (int m=0; m<nterms; m++) {
      int ms = m * (m + 1) / 2 + m;
      int mms = m * (m - 1) / 2 + m - 1;
      int mps = (m + 1) * (m + 2) / 2 + m;
      if (m > 0) Ynm[ms] = Ynm[mms] * y * Anm1[ms];
      if (m < nterms-1) Ynm[mps] = x * Ynm[ms] * Anm1[mps];
      for (int n=m+2; n<nterms; n++) {
        int nms = n * (n + 1) / 2 + m;
        int nm1 = n * (n - 1) / 2 + m;
        int nm2 = (n - 1) * (n - 2) / 2 + m;
        Ynm[nms] = Anm1[nms] * x * Ynm[nm1] - Anm2[nms] * Ynm[nm2];
      }
    }
    for (int n=0; n<nterms; n++) {
      for (int m=0; m<=n; m++) {
        int nms = n * (n + 1) / 2 + m;
        Ynm[nms] *= sqrt(2 * n + 1.0);
      }
    }
  }

  void get_Ynmd(int nterms, real_t x, real_t * Ynm, real_t * Ynmd) {
    real_t y = -sqrt((1 - x) * (1 + x));
    real_t y2 = y * y;
    Ynm[0] = 1;
    Ynmd[0] = 0;
    Ynm[1] = x * Ynm[0] * Anm1[1];
    Ynmd[1] = (x * Ynmd[0] + Ynm[0]) * Anm1[1];
    for (int n=2; n<nterms; n++) {
      int ns = n * (n + 1) / 2;
      int nm1 = n * (n - 1) / 2;
      int nm2 = (n - 1) * (n - 2) / 2;
      Ynm[ns] = Anm1[ns] * x * Ynm[nm1] - Anm2[ns] * Ynm[nm2];
      Ynmd[ns] = Anm1[ns] * (x * Ynmd[nm1] + Ynm[nm1]) - Anm2[ns] * Ynmd[nm2];
    }
    for (int m=1; m<nterms; m++) {
      int ms = m * (m + 1) / 2 + m;
      int mms = m * (m - 1) / 2 + m - 1;
      int mps = (m + 1) * (m + 2) / 2 + m;
      if (m == 1) Ynm[ms] = -Ynm[mms] * Anm1[ms];
      if (m > 1) Ynm[ms] = Ynm[mms] * y * Anm1[ms];
      if (m > 0) Ynmd[ms] = -Ynm[ms] * m * x;
      if (m < nterms-1) Ynm[mps] = x * Ynm[ms] * Anm1[mps];
      if (m < nterms-1) Ynmd[mps] = (x * Ynmd[ms] + y2 * Ynm[ms]) * Anm1[mps];
      for (int n=m+2; n<nterms; n++) {
        int nms = n * (n + 1) / 2 + m;
        int nm1 = n * (n - 1) / 2 + m;
        int nm2 = (n - 1) * (n - 2) / 2 + m;
        Ynm[nms] = Anm1[nms] * x * Ynm[nm1] - Anm2[nms] * Ynm[nm2];
        Ynmd[nms] = Anm1[nms] * (x * Ynmd[nm1] + y2 * Ynm[nm1]) - Anm2[nms] * Ynmd[nm2];
      }
    }
    for (int n=0; n<nterms; n++) {
      for (int m=0; m<=n; m++) {
        int nms = n * (n + 1) / 2 + m;
        Ynm[nms] *= sqrt(2 * n + 1.0);
        Ynmd[nms] *= sqrt(2 * n + 1.0);
      }
    }
  }

  void get_jn(int nterms, complex_t z, real_t scale, complex_t * jn, int ifder, complex_t * jnd) {
    int iscale[P+1];
    if (abs(z) < EPS) {
      jn[0] = 1;
      for (int i=1; i<nterms; i++) {
        jn[i] = 0;
      }
      if (ifder) {
        for (int i=0; i<nterms; i++) {
          jnd[i] = 0;
        }
        jnd[1] = 1.0 / (3 * scale);
      }
      return;
    }
    complex_t zinv = real_t(1.0) / z;
    jn[nterms-2] = 0;
    jn[nterms-1] = 1;
    real_t coef = 2 * nterms - 1;
    complex_t ztmp = coef * zinv;
    jn[nterms] = ztmp;
    int ntop = nterms;
    for (int i=0; i<ntop; i++) {
      iscale[i] = 0;
    }
    jn[ntop] = 0;
    jn[ntop-1] = 1;
    for (int i=ntop-1; i>0; i--) {
      coef = 2 * i + 1;
      ztmp = coef * zinv * jn[i] - jn[i+1];
      jn[i-1] = ztmp;
      if (abs(ztmp) > 1.0/EPS) {
        jn[i] *= EPS;
        jn[i-1] *= EPS;
        iscale[i] = 1;
      }
    }
    real_t scalinv = 1.0 / scale;
    coef = 1;
    for (int i=1; i<ntop; i++) {
      coef *= scalinv;
      if(iscale[i-1] == 1) coef *= EPS;
      jn[i] *= coef;
    }
    complex_t fj0 = sin(z) * zinv;
    complex_t fj1 = fj0 * zinv - cos(z) * zinv;
    if (abs(fj1) > abs(fj0)) {
      ztmp = fj1 / (jn[1] * scale);
    } else {
      ztmp = fj0 / jn[0];
    }
    for (int i=0; i<nterms; i++) {
      jn[i] *= ztmp;
    }
    if (ifder) {
      jn[nterms] *= ztmp;
      jnd[0] = -jn[1] * scale;
      for (int i=1; i<nterms; i++) {
        coef = i / (2 * i + 1.0);
        jnd[i] = coef * scalinv * jn[i-1] - (1 - coef) * scale * jn[i+1];
      }
    }
  }

  void initKernel() {
    XQUAD.resize(P);
    XQUAD2.resize(2*P);
    WQUAD.resize(P);
    WQUAD2.resize(2*P);
    NQUAD = fmax(6, P);
    NQUAD2 = fmax(6, 2*P);
    Anm1.resize((P+1)*(P+2)/2);
    Anm2.resize((P+1)*(P+2)/2);
    legendre(NQUAD, XQUAD, WQUAD);
    legendre(NQUAD2, XQUAD2, WQUAD2);
    getAnm(Anm1, Anm2);
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

  void P2M(Cell * C) {
    real_t Ynm[P*(P+1)/2];
    complex_t ephi[P], jn[P+1], jnd[P+1];
    complex_t Mnm[P*P];
    for (int n=0; n<P*P; n++) Mnm[n] = complex_t(0,0);
    real_t kscale = 2 * C->R * abs(WAVEK);
    for (Body * B=C->body; B!=C->body+C->numBodies; B++) {
      vec3 dX = B->X - C->X;
      real_t r, theta, phi;
      cart2sph(dX, r, theta, phi);
      real_t ctheta = std::cos(theta);
      ephi[1] = exp(I * phi);
      for (int n=2; n<P; n++) {
        ephi[n] = ephi[n-1] * ephi[1];
      }
      get_Ynm(P, ctheta, Ynm);
      complex_t z = WAVEK * r;
      get_jn(P, z, kscale, jn, 0, jnd);
      for (int n=0; n<P; n++) {
        jn[n] *= B->q;
      }
      for (int n=0; n<P; n++) {
        int nm = n * n + n;
        int nms = n * (n + 1) / 2;
        Mnm[nm] += Ynm[nms] * jn[n];
        for (int m=1; m<=n; m++) {
          nms = n * (n + 1) / 2 + m;
          int npm = n * n + n + m;
          int nmm = n * n + n - m;
          complex_t Ynmjn = Ynm[nms] * jn[n];
          Mnm[npm] += Ynmjn * conj(ephi[m]);
          Mnm[nmm] += Ynmjn * ephi[m];
        }
      }
    }
    for (int n=0; n<P; n++) C->M[n] += Mnm[n] * I * WAVEK;
  }
}
#endif
