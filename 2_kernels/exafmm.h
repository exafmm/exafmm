#ifndef exafmm_h
#define exafmm_h
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace exafmm {
  // Basic type definitions
  typedef double real_t;                                        //!< Floating point type is double precision
  typedef std::complex<real_t> complex_t;                       //!< Complex type

  //! Structure of bodies
  struct Body {
    real_t X[2];                                                //!< Position
    real_t q;                                                   //!< Charge
    real_t p;                                                   //!< Potential
    real_t F[2];                                                //!< Force
  };
  typedef std::vector<Body> Bodies;                             //!< Vector of bodies

  //! Structure of cells
  struct Cell {
    int NCHILD;                                                 //!< Number of child cells
    int NBODY;                                                  //!< Number of descendant bodies
    Cell * CHILD;                                               //!< Pointer to first child cell
    Body * BODY;                                                //!< Pointer to first body
    real_t X[2];                                                //!< Cell center
    real_t R;                                                   //!< Cell radius
    std::vector<complex_t> M;                                   //!< Multipole expansion coefficients
    std::vector<complex_t> L;                                   //!< Local expansion coefficients
  };
  typedef std::vector<Cell> Cells;                              //!< Vector of cells

  // Global variables
  int P;                                                        //!< Order of expansions
  int ncrit;                                                    //!< Number of bodies per leaf cell
  real_t theta;                                                 //!< Multipole acceptance criterion
  real_t dX[2];                                                 //!< Distance array

  // L2 norm of array X
  inline real_t norm(const real_t * X) {
    return X[0] * X[0] + X[1] * X[1];
  }
}
#endif
