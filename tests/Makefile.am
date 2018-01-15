include $(top_srcdir)/Makefile.am.include

KERNEL = -DEXAFMM_LAPLACE

bin_PROGRAMS = build_tree traverse
build_tree_SOURCES = build_tree.cxx
build_tree_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE
traverse_SOURCES = traverse.cxx
traverse_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE

bin_PROGRAMS += laplace_kernel laplace_ki_kernel helmholtz_kernel stokes_kernel
laplace_kernel_SOURCES = kernel.cxx
laplace_kernel_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE
laplace_ki_kernel_SOURCES = kernel.cxx
laplace_ki_kernel_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE_KI
helmholtz_kernel_SOURCES = kernel.cxx
helmholtz_kernel_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_HELMHOLTZ
stokes_kernel_SOURCES = kernel.cxx
stokes_kernel_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_STOKES

bin_PROGRAMS += laplace laplace_ki helmholtz stokes
laplace_SOURCES = fmm.cxx
laplace_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE
laplace_ki_SOURCES = fmm.cxx
laplace_ki_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE_KI
helmholtz_SOURCES = fmm.cxx
helmholtz_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_HELMHOLTZ
stokes_SOURCES = fmm.cxx
stokes_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_STOKES

TESTS = $(bin_PROGRAMS)
