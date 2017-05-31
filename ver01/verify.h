#ifndef verify_h
#define verify_h
#include <fstream>
#include "exafmm.h"

namespace exafmm {
  //! Verify results
  class Verify {
  public:
    bool verbose;                                               //!< Print to screen

    //! Get sum of scalar component of a vector of target bodies
    double getSumScalar(const Bodies & bodies) {
      double v = 0;                                             // Initialize difference
      for (size_t b=0; b<bodies.size(); b++) {                  // Loop over bodies
        v += bodies[b].p * bodies[b].q;                         //  Sum of scalar component for Laplace
      }                                                         // End loop over bodies
      return v;                                                 // Return difference
    }

    //! Get norm of scalar component of a vector of target bodies
    double getNrmScalar(const Bodies & bodies) {
      double v = 0;                                             // Initialize norm
      for (size_t b=0; b<bodies.size(); b++) {                  // Loop over bodies
        v += std::abs(bodies[b].p * bodies[b].p);               //  Norm of scalar component
      }                                                         // End loop over bodies
      return v;                                                 // Return norm
    }

    //! Get difference between scalar component of two vectors of target bodies
    double getDifScalar(const Bodies & bodies, const Bodies & bodies2) {
      double v = 0;                                             // Initialize difference
      for (size_t b=0; b<bodies.size(); b++) {                  // Loop over bodies & bodies2
        v += std::abs((bodies[b].p - bodies2[b].p) *            //  Difference of scalar component
                      (bodies[b].p - bodies2[b].p));
      }                                                         // End loop over bodies & bodies2
      return v;                                                 // Return difference
    }

    //! Get relative difference between scalar component of two vectors of target bodies
    double getRelScalar(const Bodies & bodies, const Bodies & bodies2) {
      double v = 0;                                             // Initialize difference
      for (size_t b=0; b<bodies.size(); b++) {                  // Loop over bodies & bodies2
        v += std::abs(((bodies[b].p - bodies2[b].p) * (bodies[b].p - bodies2[b].p))
		                / (bodies2[b].p * bodies2[b].p));           //  Difference of scalar component
      }                                                         // End loop over bodies & bodies2
      return v;                                                 // Return difference
    }

    //! Get norm of scalar component of a vector of target bodies
    double getNrmVector(const Bodies & bodies) {
      double v = 0;                                             // Initialize norm
      for (size_t b=0; b<bodies.size(); b++) {                  // Loop over bodies
        v += std::abs(norm(bodies[b].F));                       //  Norm of vector component
      }                                                         // End loop over bodies
      return v;                                                 // Return norm
    }

    //! Get difference between scalar component of two vectors of target bodies
    double getDifVector(const Bodies & bodies, const Bodies & bodies2) {
      double v = 0;                                             // Initialize difference
      for (size_t b=0; b<bodies.size(); b++) {                  // Loop over bodies & bodies2
        v += std::abs((bodies[b].F[0] - bodies2[b].F[0]) * (bodies[b].F[0] - bodies2[b].F[0]) +
                      (bodies[b].F[1] - bodies2[b].F[1]) * (bodies[b].F[1] - bodies2[b].F[1]));
      }                                                         // End loop over bodies & bodies2
      return v;                                                 // Return difference
    }

    //! Get difference between scalar component of two vectors of target bodies
    double getRelVector(const Bodies & bodies, const Bodies & bodies2) {
      double v = 0;                                             // Initialize difference
      for (size_t b=0; b<bodies.size(); b++) {                  // Loop over bodies & bodies2
        v += std::abs(((bodies[b].F[0] - bodies2[b].F[0]) * (bodies[b].F[0] - bodies2[b].F[0]) +
		                   (bodies[b].F[1] - bodies2[b].F[1]) * (bodies[b].F[1] - bodies2[b].F[1]))
		                  / norm(bodies2[b].F));
      }                                                         // End loop over bodies & bodies2
      return v;                                                 // Return difference
    }
  };
}
#endif
