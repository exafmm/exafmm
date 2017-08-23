#ifndef dataset_h
#define dataset_h
#include <cassert>
#include <cmath>
#include <cstdlib>
#include "exafmm.h"

namespace exafmm {
  //! Split range and return partial range
  void splitRange(int & begin, int & end, int iSplit, int numSplit) {
    assert(end > begin);
    int size = end - begin;
    int increment = size / numSplit;
    int remainder = size % numSplit;
    begin += iSplit * increment + std::min(iSplit,remainder);
    end = begin + increment;
    if (remainder > iSplit) end++;
  }

  //! Uniform distribution on [-1,1]^3 lattice
  Bodies lattice(int numBodies, int mpirank, int mpisize) {
    int nx = int(std::pow(numBodies*mpisize, 1./3));
    int ny = nx;
    int nz = nx;
    int begin = 0;
    int end = nz;
    splitRange(begin, end, mpirank, mpisize);
    int numLattice = nx * ny * (end - begin);
    Bodies bodies(numLattice);
    for (int ix=0, b=0; ix<nx; ix++) {
      for (int iy=0; iy<ny; iy++) {
        for (int iz=begin; iz<end; ++iz, ++b) {
          bodies[b].X[0] = (ix / real_t(nx-1)) * 2 - 1;
          bodies[b].X[1] = (iy / real_t(ny-1)) * 2 - 1;
          bodies[b].X[2] = (iz / real_t(nz-1)) * 2 - 1;
        }
      }
    }
    return bodies;
  }

  //! Random distribution in [-1,1]^3 cube
  Bodies cube(int numBodies, int seed, int numSplit) {
    Bodies bodies(numBodies);
    for (int i=0; i<numSplit; i++, seed++) {
      int begin = 0;
      int end = bodies.size();
      splitRange(begin, end, i, numSplit);
      srand48(seed);
      for (size_t b=begin; b!=end; ++b) {
        for (int d=0; d<3; d++) {
          bodies[b].X[d] = drand48() * 2 * M_PI - M_PI;
        }
      }
    }
    return bodies;
  }

  //! Random distribution on r = 1 sphere
  Bodies sphere(int numBodies, int seed, int numSplit) {
    Bodies bodies(numBodies);
    for (int i=0; i<numSplit; i++, seed++) {
      int begin = 0;
      int end = bodies.size();
      splitRange(begin, end, i, numSplit);
      srand48(seed);
      for (size_t b=begin; b!=end; ++b) {
        for (int d=0; d<3; d++) {
          bodies[b].X[d] = drand48() * 2 - 1;
        }
        real_t r = std::sqrt(norm(bodies[b].X));
        bodies[b].X *= M_PI / r;
      }
    }
    return bodies;
  }

  //! Random distribution on one octant of a r = 1 sphere
  Bodies octant(int numBodies, int seed, int numSplit) {
    Bodies bodies(numBodies);
    for (int i=0; i<numSplit; i++, seed++) {
      int begin = 0;
      int end = bodies.size();
      splitRange(begin, end, i, numSplit);
      srand48(seed);
      for (size_t b=begin; b!=end; ++b) {
        real_t theta = drand48() * M_PI * 0.5;
        real_t phi = drand48() * M_PI * 0.5;
        bodies[b].X[0] = 2 * M_PI * std::sin(theta) * std::cos(phi) - M_PI;
        bodies[b].X[1] = 2 * M_PI * std::sin(theta) * std::sin(phi) - M_PI;
        bodies[b].X[2] = 2 * M_PI * std::cos(theta) - M_PI;
      }
    }
    return bodies;
  }

  //! Plummer distribution in a r = M_PI/2 sphere
  Bodies plummer(int numBodies, int seed, int numSplit) {
    Bodies bodies(numBodies);
    for (int i=0; i<numSplit; i++, seed++) {
      int begin = 0;
      int end = bodies.size();
      splitRange(begin, end, i, numSplit);
      srand48(seed);
      size_t b = begin;
      while (b != end) {
        real_t X1 = drand48();
        real_t X2 = drand48();
        real_t X3 = drand48();
        real_t R = 1.0 / sqrt( (pow(X1, -2.0 / 3.0) - 1.0) );
        if (R < 100.0) {
          real_t Z = (1.0 - 2.0 * X2) * R;
          real_t X = sqrt(R * R - Z * Z) * std::cos(2.0 * M_PI * X3);
          real_t Y = sqrt(R * R - Z * Z) * std::sin(2.0 * M_PI * X3);
          real_t scale = 3.0 * M_PI / 16.0;
          X *= scale; Y *= scale; Z *= scale;
          bodies[b].X[0] = X;
          bodies[b].X[1] = Y;
          bodies[b].X[2] = Z;
          ++b;
        }
      }
    }
    return bodies;
  }

  //! Initialize source values
  void initSource(Bodies & bodies, int seed, int numSplit) {
    for (int i=0; i<numSplit; i++, seed++) {
      int begin = 0;
      int end = bodies.size();
      splitRange(begin, end, i, numSplit);
      srand48(seed);

      real_t average = 0;
      for (size_t b=begin; b!=end; ++b) {
        bodies[b].q = drand48() - .5;
        average += bodies[b].q;
      }
      average /= (end - begin);
      for (size_t b=begin; b!=end; ++b) {
        bodies[b].q -= average;
      }
    }
  }

  //! Initialize target values
  void initTarget(Bodies & bodies) {
    for (size_t b=0; b!=bodies.size(); ++b) {
      bodies[b].p = 0;
      bodies[b].F = 0;
    }
  }

  //! Initialize dsitribution, source & target value of bodies
  Bodies initBodies(int numBodies, const char * distribution,
                    int mpirank=0, int mpisize=1, int numSplit=1) {
    Bodies bodies;
    switch (distribution[0]) {
      case 'l':
        bodies = lattice(numBodies,mpirank,mpisize);
        break;
      case 'c':
        bodies = cube(numBodies,mpirank,numSplit);
        break;
      case 's':
        bodies = sphere(numBodies,mpirank,numSplit);
        break;
      case 'o':
        bodies = octant(numBodies,mpirank,numSplit);
        break;
      case 'p':
        bodies = plummer(numBodies,mpirank,numSplit);
        break;
      default:
        fprintf(stderr, "Unknown data distribution %s\n", distribution);
    }
    initSource(bodies,mpirank,numSplit);
    initTarget(bodies);
    return bodies;
  }

  //! Sample a subset of target bodies
  void sampleBodies(Bodies & bodies, int numTargets) {
    if (int(bodies.size()) > numTargets) {
      int stride = bodies.size() / numTargets;
      for (int b=0; b<numTargets; b++) {
        bodies[b] = bodies[b*stride];
      }
      bodies.resize(numTargets);
    }
  }
}
#endif
