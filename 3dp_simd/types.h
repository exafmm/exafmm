#ifndef types_h
#define types_h
#include <assert.h>                                             // Some compilers don't have cassert
#include <complex>
// #include "kahan.h"
#include "macros.h"
#include <stdint.h>
#include <vector>
#include "vec.h"

namespace exafmm {
  // Basic type definitions
  typedef vec<3,int> ivec3;                                     //!< Vector of 3 int types
  typedef vec<3,real_t> vec3;                                   //!< Vector of 3 real_t types
  typedef vec<4,real_t> vec4;                                   //!< Vector of 4 real_t types
  typedef vec<3,complex_t> cvec3;                               //!< Vector of 3 complex_t types

  // SIMD vector types for AVX512, AVX, and SSE
  const int NSIMD = SIMD_BYTES / int(sizeof(real_t));           //!< SIMD vector length (SIMD_BYTES defined in macros.h)
  typedef vec<NSIMD,real_t> simdvec;                            //!< SIMD vector type

  // Kahan summation types (Achieves quasi-double precision using single precision types)
/* #if EXAFMM_USE_KAHAN */
//   typedef kahan<real_t> kreal_t;                                //!< Real type with Kahan summation
//   typedef kahan<complex_t> kcomplex_t;                          //!< Complex type with Kahan summation
//   typedef kahan<simdvec> ksimdvec;                              //!< SIMD vector type with Kahan summation
// #else
/*   typedef real_t kreal_t;                                       //!< Real type (dummy Kahan) */
  // typedef complex_t kcomplex_t;                                 //!< Complex type (dummy Kahan)
  /* typedef simdvec ksimdvec;                                      *///!< SIMD vector type (dummy Kahan)
// #endif
/*   typedef vec<4,kreal_t> kvec4;                                 //!< Vector of 4 real types with Kahan summaiton */
  // typedef vec<4,kcomplex_t> kcvec4;                             //!< Vector of 4 complex types with Kahan summaiton

  //! Center and radius of bounding box
  struct Box {
    vec3   X;                                                   //!< Box center
    real_t R;                                                   //!< Box radius
  };

  //! Min & max bounds of bounding box
  struct Bounds {
    vec3 Xmin;                                                  //!< Minimum value of coordinates
    vec3 Xmax;                                                  //!< Maximum value of coordinates
  };

  //! Equations supported
  enum Equation {
    Empty,                                                      //!< Empty kernel
    Laplace,                                                    //!< Laplace kernel
    Helmholtz,                                                  //!< Helmholtz kernel
    BiotSavart                                                  //!< Biot-Savart kernel
  };

  //! Basis of expansion supported
  enum Basis {
    Cartesian,                                                  //! Cartesian Taylor expansion
    Spherical                                                   //! Spherical Harmonics expansion
  };

  //! Structure of aligned source for SIMD
  template<Equation equation=Empty>
  struct Source {                                               //!< Base components of source structure
    vec3      X;                                                //!< Position
  };
  template<>
  struct Source<Laplace> : public Source<> {                    //!< Special components for Laplace
    real_t    SRC;                                              //!< Scalar real values
  };
  template<>
  struct Source<Helmholtz> : public Source<> {                  //!< Special components for Helmholtz
    complex_t SRC;                                              //!< Scalar complex values
  };
  template<>
  struct Source<BiotSavart> : public Source<> {                 //!< Special components for Biot-Savart
    vec4      SRC;                                              //!< Vector real values
  };

  struct KernelBase {
    static vec3 Xperiodic;                                      //!< Periodic coordinate offset
    static real_t eps2;                                         //!< Epslion squared
    static complex_t wavek;                                     //!< Helmholtz wave number
  };

  // Multipole/local expansion coefficients
#ifdef EXAFMM_PMAX
  const int Pmax = EXAFMM_PMAX;                                 //!< Order of expansions
#else
  const int Pmax = 10;                                          //!< Order of expansions
#endif
  const int Pmin = 4;
}
#endif
