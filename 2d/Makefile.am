include $(top_srcdir)/Makefile.am.include

bin_PROGRAMS = laplace tree list 
laplace_SOURCES = fmm.cxx
laplace_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_EAGER
tree_SOURCES = tree.cxx
tree_CPPFLAGS = $(laplace_CPPFLAGS)
list_SOURCES = list.cxx
list_CPPFLAGS = $(laplace_CPPFLAGS)

run_laplace: laplace
	./$< -n 10000 -a -v
run_tree: tree
	./$< -n 100000
run_list: list
	./$< -n 100000

bin_PROGRAMS += laplace_kernel
laplace_kernel_SOURCES = kernel.cxx
laplace_kernel_CPPFLAGS = $(laplace_CPPFLAGS)

run_laplace_kernel: laplace_kernel
	./$< -P 10
	./$< -P 20
	./$< -P 30
