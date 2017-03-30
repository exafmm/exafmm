#ifndef exafmm_h
#define exafmm_h
#include <complex>
#include <cstdlib>
#include <cstdio>
#include <vector>

namespace exafmm {
  //! Basic type definitions
  typedef double real_t;                                        //!< Floating point type
  typedef std::complex<real_t> complex_t;                       //!< Complex type

  //! Structure of bodies
  struct Body {
    real_t X[3];                                                //!< Position
    real_t q;                                                   //!< Charge
    real_t p;                                                   //!< Potential
    real_t F[3];                                                //!< Force
  };
  typedef std::vector<Body> Bodies;                             //!< Vector of bodies

  //! Structure of cells
  struct Cell {
    int NCHILD;                                                 //!< Number of child cells
    int NBODY;                                                  //!< Number of descendant bodies
    Cell * CHILD;                                               //!< Pointer of first child cell
    Body * BODY;                                                //!< Pointer of first body
    real_t X[3];                                                //!< Cell center
    real_t R;                                                   //!< Cell radius
#if EXAFMM_LAZY
    std::vector<Cell*> listM2L;                                 //!< M2L interaction list
    std::vector<Cell*> listP2P;                                 //!< P2P interaction list
#endif
    std::vector<complex_t> M;                                   //!< Multipole expansion coefs
    std::vector<complex_t> L;                                   //!< Local expansion coefs
  };
  typedef std::vector<Cell> Cells;                              //!< Vector of cells

  //! Global variables
  int P;                                                        //!< Order of expansions
  int NTERM;                                                    //!< Number of coefficients
  int ncrit;                                                    //!< Number of bodies per leaf cell
  int images;                                                   //!< Number of periodic image sublevels
  real_t cycle = 2 * M_PI;                                      //!< Cycle of periodic boundary condition
  real_t theta;                                                 //!< Multipole acceptance criterion
  real_t dX[3];                                                 //!< Distance vector
  real_t Xperiodic[3];                                          //!< Periodic coordinate offset
#pragma omp threadprivate(dX,Xperiodic)                         //!< Make global variables private
}
#endif
