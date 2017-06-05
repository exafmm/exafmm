#ifndef base_mpi_h
#define base_mpi_h
#include <mpi.h>
#include "exafmm.h"

namespace exafmm {
  //! Custom MPI utilities
  class BaseMPI {
  private:
    int external;                                               //!< Flag to indicate external MPI_Init/Finalize
    const int wait;                                             //!< Waiting time between output of different ranks

  public:
    int mpirank;                                                //!< Rank of MPI communicator
    int mpisize;                                                //!< Size of MPI communicator

    //! Constructor
    BaseMPI() : external(0), wait(100) {                        // Initialize variables
      int argc(0);                                              // Dummy argument count
      char ** argv;                                             // Dummy argument value
      MPI_Initialized(&external);                               // Check if MPI_Init has been called
      if (!external) MPI_Init(&argc, &argv);                    // Initialize MPI communicator
      MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);                  // Get rank of current MPI process
      MPI_Comm_size(MPI_COMM_WORLD, &mpisize);                  // Get number of MPI processes
    }

    //! Destructor
    ~BaseMPI() {
      if (!external) MPI_Finalize();                            // Finalize MPI communicator
    }
  };

  //! Allreduce int type from all ranks
  int allreduceInt(int send) {
    int recv;                                                   // Receive buffer
    MPI_Allreduce(&send, &recv, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); // Communicate values
    return recv;                                                // Return received values
  }

  //! Allreduce bounding box from all ranks
  void allreduceBounds(real_t & R0, real_t * X0, const real_t & r0, real_t * x0 ) {
    float Xmin[2], Xmax[2];
    float xmin[2], xmax[2];
    for (int d=0; d<2; ++d) {
      xmin[d] = x0[d] - r0;
      xmax[d] = x0[d] + r0;
    }
    MPI_Allreduce(xmin, Xmin, 2, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(xmax, Xmax, 2, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    for (int d=0; d<2; ++d) X0[d] = (Xmin[d] + Xmax[d]) / 2;
    R0 = 0;
    for (int d=0; d<2; ++d) {
      R0 = fmax(X0[d] - Xmin[d], R0);
      R0 = fmax(Xmax[d] - X0[d], R0);
    }
    R0 *= 1.00001;
  }
}
#endif
