#ifndef exafmm_h
#define exafmm_h
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "vec.h"

namespace exafmm {
  //! Basic type definitions
  typedef double real_t;                                        //!< Floating point type
  typedef std::complex<real_t> complex_t;                       //!< Complex type
  typedef vec<2,real_t> vec2;                                   //!< Vector of 2 real_t types

  //! Structure of bodies
  struct Body {
    vec2 X;                                                     //!< Position
    real_t q;                                                   //!< Charge
    real_t p;                                                   //!< Potential
    vec2 F;                                                     //!< Force
  };
  typedef std::vector<Body> Bodies;                             //!< Vector of bodies

  //! Structure of cells
  struct Cell {
    int NCHILD;                                                 //!< Number of child cells
    int NBODY;                                                  //!< Number of descendant bodies
    Cell * CHILD;                                               //!< Pointer to first child cell
    Body * BODY;                                                //!< Pointer to first body
    vec2 X;                                                     //!< Cell center
    real_t R;                                                   //!< Cell radius
#if EXAFMM_LAZY
    std::vector<Cell*> listM2L;                                 //!< M2L interaction list
    std::vector<Cell*> listP2P;                                 //!< P2P interaction list
#endif
    std::vector<complex_t> M;                                   //!< Multipole expansion coefficients
    std::vector<complex_t> L;                                   //!< Local expansion coefficients
  };
  typedef std::vector<Cell> Cells;                              //!< Vector of cells

  //! Global variables
  int P;                                                        //!< Order of expansions
  int ncrit;                                                    //!< Number of bodies per leaf cell
  real_t theta;                                                 //!< Multipole acceptance criterion
}
#endif
