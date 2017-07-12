include $(top_srcdir)/Makefile.am.include

bin_PROGRAMS = build_tree traverse kernel fmm
build_tree_SOURCES = build_tree.cxx
build_tree_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_EAGER
traverse_SOURCES = traverse.cxx
traverse_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_EAGER
kernel_SOURCES = kernel.cxx
kernel_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_EAGER
fmm_SOURCES = fmm.cxx
fmm_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_EAGER

TESTS = $(bin_PROGRAMS)
