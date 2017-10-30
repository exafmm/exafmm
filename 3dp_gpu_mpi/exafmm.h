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
  MPI_Datatype MPI_REAL_T = MPI_FLOAT;          //!< Floating point MPI type is single precision
  const real_t EPS = 1e-8f;                     //!< Single precision epsilon
#else
  typedef double real_t;                        //!< Floating point type is double precision
  MPI_Datatype MPI_REAL_T = MPI_DOUBLE;         //!< Floating point MPI type is double precision
  const real_t EPS = 1e-16;                     //!< Double precision epsilon
#endif
  typedef std::complex<real_t> complex_t;       //!< Complex type
  typedef vec<3,int> ivec3;                     //!< Vector of 3 int types
  typedef vec<3,real_t> vec3;                   //!< Vector of 3 real_t types
  const complex_t I(0.,1.);                     //!< Imaginary unit

  //! Structure of bodies
  struct Body {
    vec3 X;                                     //!< Position
    real_t q;                                   //!< Charge
    real_t p;                                   //!< Potential
    vec3 F;                                     //!< Force
    uint64_t key;                               //!< Hilbert key
  };
  typedef std::vector<Body> Bodies;             //!< Vector of bodies

  //! Base components of cells
  struct CellBase {
    int numChilds;                              //!< Number of child cells
    int numBodies;                              //!< Number of descendant bodies
    vec3 X;                                     //!< Cell center
    real_t R;                                   //!< Cell radius
    uint64_t key;                               //!< Hilbert key
  };

  //! Structure of cells
  struct Cell : public CellBase {
    Cell * child;                               //!< Pointer of first child cell
    Body * body;                                //!< Pointer of first body
    std::vector<complex_t> M;                   //!< Multipole expansion coefs
    std::vector<complex_t> L;                   //!< Local expansion coefs
    using CellBase::operator=;                  //!< Substitution to derived class
  };
  typedef std::vector<Cell> Cells;              //!< Vector of cells

  //! Global variables
  int P;                                        //!< Order of expansions
  int NTERM;                                    //!< Number of coefficients
  int NCRIT;                                    //!< Number of bodies per leaf cell
  int IMAGES;                                   //!< Number of periodic image sublevels
  int IX[3];                                    //!< 3-D periodic index
  real_t CYCLE;                                 //!< Cycle of periodic boundary condition
  real_t THETA;                                 //!< Multipole acceptance criterion
  real_t R0;                                    //!< Radius of the bounding box
  vec3 X0;                                      //!< Center of the bounding box
}
#endif
