#ifndef exafmm_h
#define exafmm_h
#include <complex>
#include <cstdlib>
#include <cstdio>
#include <vector>

namespace exafmm {
  // Basic type definitions
  typedef double real_t;                                        //!< Floating point type
  typedef std::complex<real_t> complex_t;                       //!< Complex type
  const real_t EPS = 1e-16;                                     // Double precision epsilon

  //! Structure of bodies
  struct Body {
    real_t X[3];                                                //!< Position
    real_t q;                                                   //!< Charge
    real_t p;                                                   //!< Potential
    real_t F[3];                                                //!< Force
    int irank;                                                  //!< Initial rank numbering for partitioning back
    real_t weight;                                              //!< Weight for partitioning
  };
  typedef std::vector<Body> Bodies;                             //!< Vector of bodies
  typedef typename Bodies::iterator B_iter;

  //! Min & max bounds of bounding box
  struct Bounds {
    real_t Xmin[3];                                             //!< Minimum value of coordinates
    real_t Xmax[3];                                             //!< Maximum value of coordinates
  };

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
    std::vector<int> periodicM2L;                               //!< M2L periodic index
    std::vector<int> periodicP2P;                               //!< P2P periodic index
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
  int iX[3];                                                    //!< 3-D periodic index
  real_t cycle;                                                 //!< Cycle of periodic boundary condition
  real_t theta;                                                 //!< Multipole acceptance criterion
  real_t dX[3];                                                 //!< Distance vector
#pragma omp threadprivate(iX,dX)                                //!< Make global variables private

  //!< L2 norm of vector X
  inline real_t norm(const real_t * X) {
    return X[0] * X[0] + X[1] * X[1] + X[2] * X[2];             // L2 norm
  }
}
#endif
