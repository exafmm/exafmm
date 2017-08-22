#ifndef ewald_h
#define ewald_h
#include "exafmm.h"

namespace exafmm {
  //! Wave structure for Ewald summation
  struct Wave {
    vec3 K;                                     //!< 3-D wave number vector
    real_t real;                                //!< real part of wave
    real_t imag;                                //!< imaginary part of wave
  };
  typedef std::vector<Wave> Waves;              //!< Vector of Wave types

  static int KSIZE;                             //!< Number of waves in Ewald summation
  static real_t ALPHA;                          //!< Scaling parameter for Ewald summation
  static real_t SIGMA;                          //!< Scaling parameter for Ewald summation
  static real_t CUTOFF;                         //!< Cutoff distance
  static vec3 K;                                //!< Wave number vector
  static vec3 SCALE;                            //!< Scale vector

  //! Forward DFT
  void dft(Waves & waves, Bodies & bodies) {
#pragma omp parallel for
    for (size_t w=0; w<waves.size(); w++) {
      waves[w].real = waves[w].imag = 0;
      for (size_t b=0; b<bodies.size(); b++) {
        real_t th = sum(waves[w].K * bodies[b].X * SCALE);
        waves[w].real += bodies[b].q * std::cos(th);
        waves[w].imag += bodies[b].q * std::sin(th);
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
        real_t th = sum(waves[w].K * bodies[b].X * SCALE);
        real_t dtmp = waves[w].real * std::sin(th) - waves[w].imag * std::cos(th);
        p += waves[w].real * std::cos(th) + waves[w].imag * std::sin(th);
        F -= waves[w].K * dtmp;
      }
      F *= SCALE;
      bodies[b].p += p;
      bodies[b].F += F;
    }
  }

  //! Initialize wave vector
  Waves initWaves() {
    for (int d=0; d<3; d++) SCALE[d] = 2 * M_PI / CYCLE;
    Waves waves;
    int kmaxsq = KSIZE * KSIZE;
    int kmax = KSIZE;
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
            wave.real = wave.imag = 0;
            waves.push_back(wave);
          }
        }
      }
    }
    return waves;
  }

  //! Ewald real part P2P kernel
  void realP2P(Cell * Ci, Cell * Cj) {
    for (Body * Bi=Ci->body; Bi!=Ci->body+Ci->numBodies; Bi++) {
      for (Body * Bj=Cj->body; Bj!=Cj->body+Cj->numBodies; Bj++) {
        vec3 dX;
        for (int d=0; d<3; d++) dX[d] = Bi->X[d] - Bj->X[d] - IX[d] * CYCLE;
        real_t R2 = norm(dX);
        if (0 < R2 && R2 < CUTOFF * CUTOFF) {
          real_t R2s = R2 * ALPHA * ALPHA;
          real_t Rs = std::sqrt(R2s);
          real_t invRs = 1 / Rs;
          real_t invR2s = invRs * invRs;
          real_t invR3s = invR2s * invRs;
          real_t dtmp = Bj->q * (M_2_SQRTPI * std::exp(-R2s) * invR2s + erfc(Rs) * invR3s);
          dtmp *= ALPHA * ALPHA * ALPHA;
          Bi->p += Bj->q * erfc(Rs) * invRs * ALPHA;
          Bi->F -= dX * dtmp;
        }
      }
    }
  }

  //! Get leaf cells
  void getLeaf(Cells & ileafs, Cell * Ci) {
    if (Ci->numChilds == 0) ileafs.push_back(*Ci);
    for (Cell * ci=Ci->child; ci!=Ci->child+Ci->numChilds; ci++) {
      getLeaf(ileafs, ci);
    }
  }

  //! Find neighbor cell
  void neighbor(Cell * Ci, Cell * Cj) {
    vec3 dX = Ci->X - Cj->X;
    for (int d=0; d<3; d++) {
      IX[d] = 0;
      if(dX[d] < -CYCLE / 2) IX[d]--;
      if(dX[d] >  CYCLE / 2) IX[d]++;
      dX[d] -= IX[d] * CYCLE;
    }
    real_t R = std::sqrt(norm(dX));
    if (R - Ci->R - Cj->R < sqrtf(3) * CUTOFF) {
      if(Cj->numChilds == 0) realP2P(Ci, Cj);
      for (Cell * cj=Cj->child; cj!=Cj->child+Cj->numChilds; cj++) {
        neighbor(Ci, cj);
      }
    }
  }

  //! Ewald real part
  void realPart(Cells & icells, Cells & jcells) {
    Cells ileafs;
    getLeaf(ileafs, &icells[0]);
#pragma omp parallel for
    for (size_t i=0; i<ileafs.size(); i++) {
      neighbor(&ileafs[i], &jcells[0]);
    }
  }

  //! Ewald wave part
  void wavePart(Bodies & bodies, Bodies & jbodies) {
    Waves waves = initWaves();
    dft(waves,jbodies);
    real_t coef = 2 / SIGMA / CYCLE / CYCLE / CYCLE;
    real_t coef2 = 1 / (4 * ALPHA * ALPHA);
    for (size_t w=0; w<waves.size(); w++) {
      K = waves[w].K * SCALE;
      real_t K2 = K[0] * K[0] + K[1] * K[1] + K[2] * K[2];
      real_t factor = coef * std::exp(-K2 * coef2) / K2;
      waves[w].real *= factor;
      waves[w].imag *= factor;
    }
    idft(waves,bodies);
  }

  //! Subtract self term
  void selfTerm(Bodies & bodies) {
    for (size_t b=0; b<bodies.size(); b++) {
      bodies[b].p -= M_2_SQRTPI * bodies[b].q * ALPHA;
    }
  }
}
#endif
