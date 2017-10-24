include $(top_srcdir)/Makefile.am.include

KERNEL = -DEXAFMM_LAPLACE

bin_PROGRAMS = build_tree traverse
build_tree_SOURCES = build_tree.cxx
build_tree_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE
traverse_SOURCES = traverse.cxx
traverse_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE

bin_PROGRAMS += laplace_kernel helmholtz_kernel
laplace_kernel_SOURCES = kernel.cxx
laplace_kernel_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE
helmholtz_kernel_SOURCES = kernel.cxx
helmholtz_kernel_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_HELMHOLTZ

bin_PROGRAMS += laplace helmholtz
laplace_SOURCES = fmm.cxx
laplace_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE
helmholtz_SOURCES = fmm.cxx
helmholtz_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_HELMHOLTZ

TESTS = $(bin_PROGRAMS)
