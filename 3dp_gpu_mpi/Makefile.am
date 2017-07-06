include $(top_srcdir)/Makefile.am.include

bin_PROGRAMS = tree list kernel fmm
tree_SOURCES = tree.cxx
tree_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_EAGER
list_SOURCES = list.cxx
list_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_EAGER
kernel_SOURCES = kernel.cxx
kernel_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_EAGER
fmm_SOURCES = fmm.cxx
fmm_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_EAGER

TESTS = $(bin_PROGRAMS)
