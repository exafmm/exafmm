#ifndef verify_h
#define verify_h
#include <fstream>
#include "exafmm.h"

namespace exafmm {
  class Verify {
  private:
    const char * path;

  public:
    Verify(const char * _path="./") : path(_path) {}

    double getSumScalar(const Bodies & bodies) {
      double v = 0;
      for (size_t b=0; b<bodies.size(); b++) {
        v += bodies[b].p * bodies[b].q;
      }
      return v;
    }

    double getNrmScalar(const Bodies & bodies) {
      double v = 0;
      for (size_t b=0; b<bodies.size(); b++) {
        v += bodies[b].p * bodies[b].p;
      }
      return v;
    }

    double getDifScalar(const Bodies & bodies, const Bodies & bodies2) {
      double v = 0;
      for (size_t b=0; b<bodies.size(); b++) {
        v += (bodies[b].p - bodies2[b].p) * (bodies[b].p - bodies2[b].p);
      }
      return v;
    }

    double getRelScalar(const Bodies & bodies, const Bodies & bodies2) {
      double v = 0;
      for (size_t b=0; b<bodies.size(); b++) {
        v += (bodies[b].p - bodies2[b].p) * (bodies[b].p - bodies2[b].p)
          / (bodies2[b].p * bodies2[b].p);
      }
      return v;
    }

    double getNrmVector(const Bodies & bodies) {
      double v = 0;
      for (size_t b=0; b<bodies.size(); b++) {
        v += norm(bodies[b].F);
      }
      return v;
    }

    double getDifVector(const Bodies & bodies, const Bodies & bodies2) {
      double v = 0;
      for (size_t b=0; b<bodies.size(); b++) {
        v += norm(bodies[b].F - bodies2[b].F);
      }
      return v;
    }

    double getRelVector(const Bodies & bodies, const Bodies & bodies2) {
      double v = 0;
      for (size_t b=0; b<bodies.size(); b++) {
        v += norm(bodies[b].F - bodies2[b].F) / norm(bodies2[b].F);
      }
      return v;
    }
  };
}
#endif
