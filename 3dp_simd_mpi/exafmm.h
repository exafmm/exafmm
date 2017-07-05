#ifndef exafmm_h
#define exafmm_h
#include <complex>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include "vec.h"

namespace exafmm {
  // Basic type definitions
  typedef double real_t;                                        //!< Floating point type is double precision
  typedef std::complex<real_t> complex_t;                       //!< Complex type
  typedef vec<3,real_t> vec3;                                   //!< Vector of 3 real_t types

  //! Structure of bodies
  struct Body {
    vec3 X;                                                     //!< Position
    real_t q;                                                   //!< Charge
    real_t p;                                                   //!< Potential
    vec3 F;                                                     //!< Force
  };
  typedef std::vector<Body> Bodies;                             //!< Vector of bodies

  //! Structure of cells
  struct Cell {
    int NCHILD;                                                 //!< Number of child cells
    int NBODY;                                                  //!< Number of descendant bodies
    Cell * CHILD;                                               //!< Pointer of first child cell
    Body * BODY;                                                //!< Pointer of first body
    vec3 X;                                                     //!< Cell center
    real_t R;                                                   //!< Cell radius
#if EXAFMM_LAZY
    std::vector<Cell*> listM2L;                                 //!< M2L interaction list
    std::vector<Cell*> listP2P;                                 //!< P2P interaction list
    std::vector<int> periodicM2L;                               //!< M2L periodic index
    std::vector<int> periodicP2P;                               //!< P2P periodic index
#endif
    std::vector<complex_t> M;                                   //!< Multipole expansion coefs
    std::vector<complex_t> L;                                   //!< Local expansion coefs
  };
  typedef std::vector<Cell> Cells;                              //!< Vector of cells

  // Global variables
  int P;                                                        //!< Order of expansions
  int NTERM;                                                    //!< Number of coefficients
  int ncrit;                                                    //!< Number of bodies per leaf cell
  int images;                                                   //!< Number of periodic image sublevels
  int iX[3];                                                    //!< 3-D periodic index
  real_t cycle;                                                 //!< Cycle of periodic boundary condition
  real_t theta;                                                 //!< Multipole acceptance criterion
#pragma omp threadprivate(iX)                                   //!< Make global variables private
}
#endif
