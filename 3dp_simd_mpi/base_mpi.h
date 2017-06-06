#ifndef base_mpi_h
#define base_mpi_h
#include <mpi.h>
#include "exafmm.h"
#include <iostream>
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

  //! Allreduce bounds type from all ranks
  Bounds allreduceBounds(Bounds local) {
    float localXmin[3], localXmax[3], globalXmin[3], globalXmax[3];
    for (int d=0; d<3; d++) {                                 // Loop over dimensions
      localXmin[d] = local.Xmin[d];                           //  Convert Xmin to float
      localXmax[d] = local.Xmax[d];                           //  Convert Xmax to float
    }                                                         // End loop over dimensions
    MPI_Allreduce(localXmin, globalXmin, 3, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);// Reduce domain Xmin
    MPI_Allreduce(localXmax, globalXmax, 3, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);// Reduce domain Xmax
    Bounds global;
    for (int d=0; d<3; d++) {                                 // Loop over dimensions
      real_t leeway = (globalXmax[d] - globalXmin[d]) * 1e-6; //  Adding a bit of leeway to global domain
      global.Xmin[d] = globalXmin[d] - leeway;                //  Convert Xmin to real_t
      global.Xmax[d] = globalXmax[d] + leeway;                //  Convert Xmax to real_t
    }                                                         // End loop over dimensions
    return global;                                            // Return global bounds
  }
}
#endif
