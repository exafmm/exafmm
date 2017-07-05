#ifndef ewald_h
#define ewald_h
#include "exafmm.h"

namespace exafmm {
  //! Wave structure for Ewald summation
  struct Wave {
    real_t K[3];                                                //!< 3-D wave number vector
    real_t REAL;                                                //!< real part of wave
    real_t IMAG;                                                //!< imaginary part of wave
  };
  typedef std::vector<Wave> Waves;                              //!< Vector of Wave types

  static int ksize;                                             //!< Number of waves in Ewald summation
  static real_t alpha;                                          //!< Scaling parameter for Ewald summation
  static real_t sigma;                                          //!< Scaling parameter for Ewald summation
  static real_t cutoff;                                         //!< Cutoff distance
  static real_t K[3];                                           //!< Wave number vector
  static real_t scale[3];                                       //!< Scale vector

  //! Forward DFT
  void dft(Waves & waves, Bodies & bodies) {
#pragma omp parallel for
    for (size_t w=0; w<waves.size(); w++) {
      waves[w].REAL = waves[w].IMAG = 0;
      for (size_t b=0; b<bodies.size(); b++) {
        real_t th = 0;
        for (int d=0; d<3; d++) th += waves[w].K[d] * bodies[b].X[d] * scale[d];
        waves[w].REAL += bodies[b].q * std::cos(th);
        waves[w].IMAG += bodies[b].q * std::sin(th);
      }
    }
  }

  //! Inverse DFT
  void idft(Waves & waves, Bodies & bodies) {
#pragma omp parallel for
    for (size_t b=0; b<bodies.size(); b++) {
      real_t p = 0;
      vec3 F = 0;
      for (size_t w=0; w<waves.size(); w++) {
        real_t th = 0;
        for (int d=0; d<3; d++) th += waves[w].K[d] * bodies[b].X[d] * scale[d];
        real_t dtmp = waves[w].REAL * std::sin(th) - waves[w].IMAG * std::cos(th);
        p += waves[w].REAL * std::cos(th) + waves[w].IMAG * std::sin(th);
        for (int d=0; d<3; d++) F[d] -= dtmp * waves[w].K[d];
      }
      for (int d=0; d<3; d++) F[d] *= scale[d];
      bodies[b].p += p;
      bodies[b].F += F;
    }
  }

  //! Initialize wave vector
  Waves initWaves() {
    for (int d=0; d<3; d++) scale[d]= 2 * M_PI / cycle;
    Waves waves;
    int kmaxsq = ksize * ksize;
    int kmax = ksize;
    for (int l=0; l<=kmax; l++) {
      int mmin = -kmax;
      if (l==0) mmin = 0;
      for (int m=mmin; m<=kmax; m++) {
        int nmin = -kmax;
        if (l==0 && m==0) nmin=1;
        for (int n=nmin; n<=kmax; n++) {
          real_t ksq = l * l + m * m + n * n;
          if (ksq <= kmaxsq) {
            Wave wave;
            wave.K[0] = l;
            wave.K[1] = m;
            wave.K[2] = n;
            wave.REAL = wave.IMAG = 0;
            waves.push_back(wave);
          }
        }
      }
    }
    return waves;
  }

  //! Ewald real part P2P kernel
  void realP2P(Cell * Ci, Cell * Cj) {
    for (Body * Bi=Ci->BODY; Bi!=Ci->BODY+Ci->NBODY; Bi++) {
      for (Body * Bj=Cj->BODY; Bj!=Cj->BODY+Cj->NBODY; Bj++) {
        vec3 dX;
        for (int d=0; d<3; d++) dX[d] = Bi->X[d] - Bj->X[d] - iX[d] * cycle;
        real_t R2 = norm(dX);
        if (0 < R2 && R2 < cutoff * cutoff) {
          real_t R2s = R2 * alpha * alpha;
          real_t Rs = std::sqrt(R2s);
          real_t invRs = 1 / Rs;
          real_t invR2s = invRs * invRs;
          real_t invR3s = invR2s * invRs;
          real_t dtmp = Bj->q * (M_2_SQRTPI * std::exp(-R2s) * invR2s + erfc(Rs) * invR3s);
          dtmp *= alpha * alpha * alpha;
          Bi->p += Bj->q * erfc(Rs) * invRs * alpha;
          Bi->F -= dX * dtmp;
        }
      }
    }
  }

  //! Find neighbor cell
  void neighbor(Cell * Ci, Cell * Cj) {
    vec3 dX = Ci->X - Cj->X;
    for (int d=0; d<3; d++) {
      iX[d] = 0;
      if(dX[d] < -cycle / 2) iX[d]--;
      if(dX[d] >  cycle / 2) iX[d]++;
      dX[d] -= iX[d] * cycle;
    }
    real_t R = std::sqrt(norm(dX));
    if (R - Ci->R - Cj->R < sqrtf(3) * cutoff) {
      if(Cj->NCHILD == 0) realP2P(Ci, Cj);
      for (Cell * cj=Cj->CHILD; cj!=Cj->CHILD+Cj->NCHILD; cj++) {
        neighbor(Ci, cj);
      }
    }
  }

  //! Ewald real part
  void realPart(Cell * Ci, Cell * Cj) {
    if (Ci->NCHILD == 0) neighbor(Ci, Cj);
    for (Cell * ci=Ci->CHILD; ci!=Ci->CHILD+Ci->NCHILD; ci++) {
      realPart(ci, Cj);
    }
  }

  //! Ewald wave part
  void wavePart(Bodies & bodies, Bodies & jbodies) {
    Waves waves = initWaves();
    dft(waves,jbodies);
    real_t coef = 2 / sigma / cycle / cycle / cycle;
    real_t coef2 = 1 / (4 * alpha * alpha);
    for (size_t w=0; w<waves.size(); w++) {
      for (int d=0; d<3; d++) K[d] = waves[w].K[d] * scale[d];
      real_t K2 = K[0] * K[0] + K[1] * K[1] + K[2] * K[2];
      real_t factor = coef * std::exp(-K2 * coef2) / K2;
      waves[w].REAL *= factor;
      waves[w].IMAG *= factor;
    }
    idft(waves,bodies);
  }

  //! Subtract self term
  void selfTerm(Bodies & bodies) {
    for (size_t b=0; b<bodies.size(); b++) {
      bodies[b].p -= M_2_SQRTPI * bodies[b].q * alpha;
    }
  }
}
#endif
