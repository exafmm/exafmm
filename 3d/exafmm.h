#ifndef exafmm_h
#define exafmm_h
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstdio>
#include <stdint.h>
#include <vector>
#include "vec.h"

namespace exafmm {
  //! Basic type definitions
#if EXAFMM_SINGLE
  typedef float real_t;                         //!< Floating point type is single precision
  const real_t EPS = 1e-8f;                     //!< Single precision epsilon
#else
  typedef double real_t;                        //!< Floating point type is double precision
  const real_t EPS = 1e-16;                     //!< Double precision epsilon
#endif
  typedef std::complex<real_t> complex_t;       //!< Complex type
  typedef vec<3,real_t> vec3;                   //!< Vector of 3 real_t types
#if EXAFMM_HELMHOLTZ
  typedef vec<3,complex_t> cvec3;               //!< Vector of 3 complex_t types
  const complex_t I(0.,1.);                     //!< Imaginary unit
#endif

  //! Structure of bodies
  struct Body {
    vec3 X;                                     //!< Position
#if EXAFMM_LAPLACE
    real_t q;                                   //!< Charge
    real_t p;                                   //!< Potential
    vec3 F;                                     //!< Force
#elif EXAFMM_HELMHOLTZ
    complex_t q;                                //!< Charge
    complex_t p;                                //!< Potential
    cvec3 F;                                    //!< Force
#endif
  };
  typedef std::vector<Body> Bodies;             //!< Vector of bodies

  //! Structure of cells
  struct Cell {
    int numChilds;                              //!< Number of child cells
    int numBodies;                              //!< Number of descendant bodies
    Cell * child;                               //!< Pointer of first child cell
    Body * body;                                //!< Pointer of first body
    vec3 X;                                     //!< Cell center
    real_t R;                                   //!< Cell radius
#if EXAFMM_LAZY
    std::vector<Cell*> listM2L;                 //!< M2L interaction list
    std::vector<Cell*> listP2P;                 //!< P2P interaction list
#endif
    std::vector<complex_t> M;                   //!< Multipole expansion coefs
    std::vector<complex_t> L;                   //!< Local expansion coefs
  };
  typedef std::vector<Cell> Cells;              //!< Vector of cells

  //! Global variables
  int P;                                        //!< Order of expansions
  int NTERM;                                    //!< Number of coefficients
  int NCRIT;                                    //!< Number of bodies per leaf cell
  real_t THETA;                                 //!< Multipole acceptance criterion
  complex_t WAVEK;                              //!< Wave number
}
#endif
