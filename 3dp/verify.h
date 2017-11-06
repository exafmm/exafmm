#ifndef verify_h
#define verify_h
#include "exafmm.h"

namespace exafmm {
#if !EXAFMM_STOKES
  double getSumScalar(const Bodies & bodies) {
    double v = 0;
    for (size_t b=0; b<bodies.size(); b++) {
#if EXAFMM_LAPLACE || EXAFMM_LAPLACE_KI
      v += bodies[b].p * bodies[b].q;
#elif EXAFMM_HELMHOLTZ
      v += std::abs(bodies[b].p * bodies[b].q);
#endif
    }
    return v;
  }

  double getNrmScalar(const Bodies & bodies) {
    double v = 0;
    for (size_t b=0; b<bodies.size(); b++) {
      v += std::abs(bodies[b].p * bodies[b].p);
    }
    return v;
  }

  double getDifScalar(const Bodies & bodies, const Bodies & bodies2) {
    double v = 0;
    for (size_t b=0; b<bodies.size(); b++) {
      v += std::abs((bodies[b].p - bodies2[b].p) * (bodies[b].p - bodies2[b].p));
    }
    return v;
  }
#endif

  double getNrmVector(const Bodies & bodies) {
    double v = 0;
    for (size_t b=0; b<bodies.size(); b++) {
#if EXAFMM_LAPLACE || EXAFMM_LAPLACE_KI || EXAFMM_STOKES
      v += norm(bodies[b].F);
#elif EXAFMM_HELMHOLTZ
      for (int d=0; d<3; d++) v+= std::norm(bodies[b].F[d]);
#endif
    }
    return v;
  }

  double getDifVector(const Bodies & bodies, const Bodies & bodies2) {
    double v = 0;
    for (size_t b=0; b<bodies.size(); b++) {
#if EXAFMM_LAPLACE || EXAFMM_LAPLACE_KI || EXAFMM_STOKES
      v += norm(bodies[b].F - bodies2[b].F);
#elif EXAFMM_HELMHOLTZ
      for (int d=0; d<3; d++) v+= std::norm(bodies[b].F[d] - bodies2[b].F[d]);
#endif
    }
    return v;
  }
}
#endif
