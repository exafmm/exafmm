include $(top_srcdir)/Makefile.am.include

KERNEL = -DEXAFMM_LAPLACE

bin_PROGRAMS = build_tree traverse kernel fmm
build_tree_SOURCES = build_tree.cxx
build_tree_CPPFLAGS = $(AM_CPPFLAGS) $(KERNEL)
traverse_SOURCES = traverse.cxx
traverse_CPPFLAGS = $(AM_CPPFLAGS) $(KERNEL)
kernel_SOURCES = kernel.cxx
kernel_CPPFLAGS = $(AM_CPPFLAGS) $(KERNEL)
fmm_SOURCES = fmm.cxx
fmm_CPPFLAGS = $(AM_CPPFLAGS) $(KERNEL)

TESTS = $(bin_PROGRAMS)
