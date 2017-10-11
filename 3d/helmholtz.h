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

  void rotate(real_t theta, int nterms, const complex_t * YnmIn, complex_t * YnmOut) {
    real_t Rnm1[P][2*P];
    real_t Rnm2[P][2*P];
    real_t sqrtCnm[2*P][2];
    for (int m=0; m<2*P; m++) {
      sqrtCnm[m][0] = sqrt(m+0.0);
    }
    sqrtCnm[0][1] = 0;
    sqrtCnm[1][1] = 0;
    for (int m=2; m<2*P; m++) {
      sqrtCnm[m][1] = sqrt(m * (m - 1) / 2.0);
    }
    real_t ctheta = std::cos(theta);
    if (fabs(ctheta) < EPS) ctheta = 0;
    real_t stheta = std::sin(-theta);
    if (fabs(stheta) < EPS) stheta = 0;
    real_t hsthta = stheta / sqrt(2.0);
    real_t cthtap = sqrt(2.0) * std::cos(theta * .5) * std::cos(theta * .5);
    real_t cthtan =-sqrt(2.0) * std::sin(theta * .5) * std::sin(theta * .5);
    Rnm1[0][P] = 1;
    YnmOut[0] = YnmIn[0];
    for (int n=1; n<nterms; n++) {
      for (int m=-n; m<0; m++) {
        Rnm2[0][P+m] = -sqrtCnm[n-m][1] * Rnm1[0][P+m+1];
        if (m > (1 - n)) {
          Rnm2[0][P+m] += sqrtCnm[n+m][1] * Rnm1[0][P+m-1];
        }
        Rnm2[0][P+m] *= hsthta;
        if (m > -n) {
          Rnm2[0][P+m] += Rnm1[0][P+m] * ctheta * sqrtCnm[n+m][0] * sqrtCnm[n-m][0];
        }
        Rnm2[0][P+m] /= n;
      }
      Rnm2[0][P] = Rnm1[0][P] * ctheta;
      if (n > 1) {
        Rnm2[0][P] += hsthta * sqrtCnm[n][1] * (2 * Rnm1[0][P-1]) / n;
      }
      for (int m=1; m<=n; m++) {
        Rnm2[0][P+m] = Rnm2[0][P-m];
        if (m % 2 == 0) {
          Rnm2[m][P] = Rnm2[0][P+m];
        } else {
          Rnm2[m][P] =-Rnm2[0][P+m];
        }
      }
      for (int mp=1; mp<=n; mp++) {
        real_t scale = 1 / (sqrt(2.0) * sqrtCnm[n+mp][1]);
        for (int m=mp; m<=n; m++) {
          Rnm2[mp][P+m] = Rnm1[mp-1][P+m-1] * cthtap * sqrtCnm[n+m][1];
          Rnm2[mp][P-m] = Rnm1[mp-1][P-m+1] * cthtan * sqrtCnm[n+m][1];
          if (m < (n - 1)) {
            Rnm2[mp][P+m] -= Rnm1[mp-1][P+m+1] * cthtan * sqrtCnm[n-m][1];
            Rnm2[mp][P-m] -= Rnm1[mp-1][P-m-1] * cthtap * sqrtCnm[n-m][1];
          }
          if (m < n) {
            real_t d = stheta * sqrtCnm[n+m][0] * sqrtCnm[n-m][0];
            Rnm2[mp][P+m] += Rnm1[mp-1][P+m] * d;
            Rnm2[mp][P-m] += Rnm1[mp-1][P-m] * d;
          }
          Rnm2[mp][P+m] *= scale;
          Rnm2[mp][P-m] *= scale;
          if (m > mp) {
            if ((mp+m) % 2 == 0) {
              Rnm2[m][P+mp] = Rnm2[mp][P+m];
              Rnm2[m][P-mp] = Rnm2[mp][P-m];
            } else {
              Rnm2[m][P+mp] =-Rnm2[mp][P+m];
              Rnm2[m][P-mp] =-Rnm2[mp][P-m];
            }
          }
        }
      }
      for (int m=-n; m<=n; m++) {
        int nn = n * n + n;
        int nm = n * n + n + m;
        YnmOut[nm] = YnmIn[nn] * Rnm2[0][P+m];
        for (int mp=1; mp<=n; mp++) {
          nm = n * n + n + m;
          int npm = n * n + n + mp;
          int nmm = n * n + n - mp;
          YnmOut[nm] += YnmIn[npm] * Rnm2[mp][P+m] + YnmIn[nmm] * Rnm2[mp][P-m];
        }
      }
      for (int m=-n; m<=n; m++) {
        for (int mp=0; mp<=n; mp++) {
          Rnm1[mp][P+m] = Rnm2[mp][P+m];
        }
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

  void get_hn(int nterms, complex_t z, real_t scale, complex_t * hn) {
    if (abs(z) < EPS) {
      for (int i=0; i<nterms; i++) {
        hn[i] = 0;
      }
      return;
    }
    complex_t zi = I * z;
    complex_t zinv = scale / z;
    hn[0] = exp(zi) / zi;
    hn[1] = hn[0] * (zinv - I * scale);
    real_t scale2 = scale * scale;
    for (int i=2; i<nterms; i++) {
      hn[i] = zinv * real_t(2 * i - 1.0) * hn[i-1] - scale2 * hn[i-2];
    }
  }

  void get_hnd(int nterms, complex_t z, real_t scale, complex_t * hn, complex_t * hnd) {
    if (abs(z) < EPS) {
      for (int i=0; i<nterms; i++) {
        hn[i] = 0;
        hnd[i] = 0;
      }
      return;
    }
    complex_t zi = I * z;
    complex_t zinv = real_t(1.0) / z;
    hn[0] = exp(zi) / zi;
    hn[1] = hn[0] * (zinv - I) * scale;
    hnd[0] = -hn[1] / scale;
    hnd[1] = -zinv * real_t(2.0) * hn[1] + scale * hn[0];
    for (int i=2; i<nterms; i++) {
      hn[i] = (zinv * real_t(2 * i - 1.0) * hn[i-1] - scale * hn[i-2]) * scale;
      hnd[i] = -zinv * real_t(i + 1.0) * hn[i] + scale * hn[i-1];
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

  void M2M(Cell * Ci) {
    real_t Ynm[P*(P+1)/2];
    complex_t phitemp[2*P], hn[P], ephi[2*P];
    complex_t Mnm[P*P], Mrot[P*P];
    for (int n=0; n<P*P; n++) Mnm[n] = Mrot[n] = complex_t(0,0);
    real_t kscalei = 2 * Ci->R * abs(WAVEK);
    for (Cell * Cj=Ci->child; Cj!=Ci->child+Ci->numChilds; Cj++) {
      real_t kscalej = 2 * Cj->R * abs(WAVEK);
      real_t radius = 2 * Cj->R * sqrt(3.0);
      vec3 dX = Ci->X - Cj->X;
      real_t r, theta, phi;
      cart2sph(dX, r, theta, phi);
      ephi[P+1] = exp(I * phi);
      ephi[P] = 1;
      ephi[P-1] = conj(ephi[P+1]);
      for (int n=2; n<P; n++) {
        ephi[P+n] = ephi[P+n-1] * ephi[P+1];
        ephi[P-n] = conj(ephi[P+n]);
      }
      for (int n=0; n<P; n++) {
        for (int m=-n; m<=n; m++) {
          int nm = n * n + n + m;
          Mnm[nm] = Cj->M[nm] * ephi[P+m];
        }
      }
      rotate(theta, P, Mnm, Mrot);
      for (int n=0; n<P; n++) {
        for (int m=-n; m<=n; m++) {
          int nm = n * n + n + m;
          Mnm[nm] = 0;
        }
      }
      for (int l=0; l<NQUAD2; l++) {
        real_t ctheta = XQUAD2[l];
        real_t stheta = sqrt(1 - ctheta * ctheta);
        real_t rj = (r + radius * ctheta) * (r + radius * ctheta) + (radius * stheta) * (radius * stheta);
        rj = sqrt(rj);
        real_t cthetaj = (r + radius * ctheta) / rj;
        complex_t z = WAVEK * rj;
        get_Ynm(P, cthetaj, Ynm);
        get_hn(P, z, kscalej, hn);
        for (int m=-P+1; m<P; m++) {
          int mabs = abs(m);
          phitemp[P+m] = 0;
          for (int n=mabs; n<P; n++) {
            int nm = n * n + n + m;
            int nms = n * (n + 1) / 2 + mabs;
            phitemp[P+m] += Mrot[nm] * hn[n] * Ynm[nms];
          }
        }
        get_Ynm(P, XQUAD2[l], Ynm);
        for (int m=-P+1; m<P; m++) {
          int mabs = abs(m);
          z = phitemp[P+m] * WQUAD2[l] * real_t(.5);
          for (int n=mabs; n<P; n++) {
            int nm = n * n + n + m;
            int nms = n * (n + 1) / 2 + mabs;
            Mnm[nm] += z * Ynm[nms];
          }
        }
      }
      complex_t z = WAVEK * radius;
      get_hn(P, z, kscalei, hn);
      for (int n=0; n<P; n++) {
        for (int m=-n; m<=n; m++) {
          int nm = n * n + n + m;
          Mnm[nm] /= hn[n];
        }
      }
      rotate(-theta, P, Mnm, Mrot);
      for (int n=0; n<P; n++) {
        for (int m=-n; m<=n; m++) {
          int nm = n * n + n + m;
          Mnm[nm] = ephi[P-m] * Mrot[nm];
        }
      }
      for (int n=0; n<P*P; n++) Ci->M[n] += Mnm[n];
    }
  }

  void M2L(Cell * Ci, Cell * Cj) {
    real_t Ynm[P*(P+1)/2], Ynmd[P*(P+1)/2];
    complex_t phitemp[2*P], phitempn[2*P];
    complex_t hn[P], hnd[P], jn[P+1], jnd[P+1], ephi[2*P];
    complex_t Mnm[P*P], Mrot[P*P], Lnm[P*P], Lrot[P*P], Lnmd[P*P];
    for (int n=0; n<P*P; n++) Lnm[n] = Lrot[n] = complex_t(0,0);
    real_t kscalej = 2 * Cj->R * abs(WAVEK);
    real_t kscalei = 2 * Ci->R * abs(WAVEK);
    real_t radius = 2 * Cj->R * sqrt(3.0) * .5;
    vec3 dX = Ci->X - Cj->X;
    real_t r, theta, phi;
    cart2sph(dX, r, theta, phi);
    dX /= 2 * Cj->R;
    if (fabs(dX[0]) > EPS) dX[0] = fabs(dX[0]) - .5;
    if (fabs(dX[1]) > EPS) dX[1] = fabs(dX[1]) - .5;
    if (fabs(dX[2]) > EPS) dX[2] = fabs(dX[2]) - .5;
    real_t rr = sqrt(norm(dX));
    real_t coef1 = P * 1.65 - 15.5;
    real_t coef2 = P * 0.25 + 3.0;
    int Popt = coef1 / (rr * rr) + coef2;
    assert(0 < Popt);
    assert(Popt <= 2*P);
    if(Popt > P) Popt = P;
    ephi[P+1] = exp(I * phi);
    ephi[P] = 1;
    ephi[P-1] = conj(ephi[P+1]);
    for (int n=2; n<P; n++) {
      ephi[P+n] = ephi[P+n-1] * ephi[P+1];
      ephi[P-n] = conj(ephi[P+n]);
    }
    for (int n=0; n<Popt; n++) {
      for (int m=-n; m<=n; m++) {
        int nm = n * n + n + m;
        Mnm[nm] = Cj->M[nm] * ephi[P+m];
      }
    }
    rotate(theta, Popt, Mnm, Mrot);
    for (int l=0; l<NQUAD; l++) {
      real_t ctheta = XQUAD[l];
      real_t stheta = sqrt(1 - ctheta * ctheta);
      real_t rj = (r + radius * ctheta) * (r + radius * ctheta) + (radius * stheta) * (radius * stheta);
      rj = sqrt(rj);
      real_t cthetaj = (r + radius * ctheta) / rj;
      real_t sthetaj = sqrt(1 - cthetaj * cthetaj);
      real_t rn = sthetaj * stheta + cthetaj * ctheta;
      real_t thetan = (cthetaj * stheta - ctheta * sthetaj) / rj;
      complex_t z = WAVEK * rj;
      get_Ynmd(Popt, cthetaj, Ynm, Ynmd);
      get_hnd(Popt, z, kscalej, hn, hnd);
      for (int n=0; n<Popt; n++) {
        hnd[n] *= WAVEK;
      }
      for (int n=1; n<Popt; n++) {
        for (int m=1; m<=n; m++) {
          int nms = n * (n + 1) / 2 + m;
          Ynm[nms] *= sthetaj;
        }
      }
      for (int m=-Popt+1; m<Popt; m++) {
        phitemp[Popt+m] = 0;
        phitempn[Popt+m] = 0;
      }
      phitemp[Popt] = Mrot[0] * hn[0];
      phitempn[Popt] = Mrot[0] * hnd[0] * rn;
      for (int n=1; n<Popt; n++) {
        int nm = n * n + n;
        int nms = n * (n + 1) / 2;
        phitemp[Popt] += Mrot[nm] * hn[n] * Ynm[nms];
        complex_t ut1 = hnd[n] * rn;
        complex_t ut2 = hn[n] * thetan;
        complex_t ut3 = ut1 * Ynm[nms] - ut2 * Ynmd[nms] * sthetaj;
        phitempn[Popt] += ut3 * Mrot[nm];
        for (int m=1; m<=n; m++) {
          nms = n * (n + 1) / 2 + m;
          int npm = n * n + n + m;
          int nmm = n * n + n - m;
          z = hn[n] * Ynm[nms];
          phitemp[Popt+m] += Mrot[npm] * z;
          phitemp[Popt-m] += Mrot[nmm] * z;
          ut3 = ut1 * Ynm[nms] - ut2 * Ynmd[nms];
          phitempn[Popt+m] += ut3 * Mrot[npm];
          phitempn[Popt-m] += ut3 * Mrot[nmm];
        }
      }
      get_Ynm(Popt, XQUAD[l], Ynm);
      for (int m=-Popt+1; m<Popt; m++) {
        int mabs = abs(m);
        z = phitemp[Popt+m] * WQUAD[l] * real_t(.5);
        for (int n=mabs; n<Popt; n++) {
          int nm = n * n + n + m;
          int nms = n * (n + 1) / 2 + mabs;
          Lnm[nm] += z * Ynm[nms];
        }
        z = phitempn[Popt+m] * WQUAD[l] * real_t(.5);
        for (int n=mabs; n<Popt; n++) {
          int nm = n * n + n + m;
          int nms = n * (n + 1) / 2 + mabs;
          Lnmd[nm] += z * Ynm[nms];
        }
      }
    }
    complex_t z = WAVEK * radius;
    get_jn(Popt, z, kscalei, jn, 1, jnd);
    for (int n=0; n<Popt; n++) {
      for (int m=-n; m<=n; m++) {
        int nm = n * n + n + m;
        complex_t zh = jn[n];
        complex_t zhn = jnd[n] * WAVEK;
        z = zh * zh + zhn * zhn;
        Lnm[nm] = (zh * Lnm[nm] + zhn * Lnmd[nm]) / z;
      }
    }
    rotate(-theta, Popt, Lnm, Lrot);
    for (int n=0; n<Popt; n++) {
      for (int m=-n; m<=n; m++) {
        int nm = n * n + n + m;
        Lnm[nm] = ephi[P-m] * Lrot[nm];
      }
    }
    for (int n=0; n<P*P; n++) Ci->L[n] += Lnm[n];
  }
}
#endif
