include ../Makefile.am.include

bin_PROGRAMS = laplace helmholtz biotsavart tree vec
laplace_SOURCES = fmm.cxx
laplace_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE -DEXAFMM_PMAX=10
helmholtz_SOURCES = fmm.cxx
helmholtz_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_HELMHOLTZ -DEXAFMM_PMAX=20
biotsavart_SOURCES = fmm.cxx
biotsavart_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_BIOTSAVART -DEXAFMM_PMAX=10
tree_SOURCES = tree.cxx
tree_CPPFLAGS = $(laplace_CPPFLAGS)
vec_SOURCES = vec.cxx
vec_CPPFLAGS = $(laplace_CPPFLAGS)

run_laplace: laplace
	./$< -aDgv -n 100000 -r 1
run_helmholtz: helmholtz
	./$< -aDgv -n 10000 -r 1 -P 20
run_biotsavart: biotsavart
	./$< -agv -n 10000 -r 1
run_tree: tree
	./$< -av -n 1000000 -r 1
run_vec: vec
	./$<

bin_PROGRAMS += laplace_kernel helmholtz_kernel biotsavart_kernel
laplace_kernel_SOURCES = kernel.cxx
laplace_kernel_CPPFLAGS = $(laplace_CPPFLAGS)
helmholtz_kernel_SOURCES = kernel.cxx
helmholtz_kernel_CPPFLAGS = $(helmholtz_CPPFLAGS)
biotsavart_kernel_SOURCES = kernel.cxx
biotsavart_kernel_CPPFLAGS = $(biotsavart_CPPFLAGS)

run_laplace_kernel: laplace_kernel
	./$< -v -P 10
	./$< -v -P 20
	./$< -v -P 30
run_helmholtz_kernel: helmholtz_kernel
	./$< -v -P 10
	./$< -v -P 20
	./$< -v -P 30
run_biotsavart_kernel: biotsavart_kernel
	./$< -v -P 10
	./$< -v -P 20
	./$< -v -P 30

if EXAFMM_HAVE_MPI
bin_PROGRAMS += laplace_mpi helmholtz_mpi biotsavart_mpi ewald_mpi
laplace_mpi_SOURCES = fmm_mpi.cxx
laplace_mpi_CPPFLAGS = $(laplace_CPPFLAGS)
helmholtz_mpi_SOURCES = fmm_mpi.cxx
helmholtz_mpi_CPPFLAGS = $(helmholtz_CPPFLAGS)
biotsavart_mpi_SOURCES = fmm_mpi.cxx
biotsavart_mpi_CPPFLAGS = $(biotsavart_CPPFLAGS)
ewald_mpi_SOURCES = ewald.cxx
ewald_mpi_CPPFLAGS = $(laplace_CPPFLAGS)

run_laplace_mpi: laplace_mpi
	$(MPIRUN) -n 2 ./$< -aDgv -n 100000 -r 1
run_helmholtz_mpi: helmholtz_mpi
	$(MPIRUN) -n 2 ./$< -aDgv -n 10000 -r 1 -P 20
run_biotsavart_mpi: biotsavart_mpi
	$(MPIRUN) -n 2 ./$< -agv -n 10000 -r 1
run_ewald_mpi: ewald_mpi
	$(MPIRUN) -n 2 ./$< -aDgv -r 1
endif
