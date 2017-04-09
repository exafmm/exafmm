# https://www.gnu.org/software/autoconf-archive/The-Macros.html#The-Macros
AC_INIT([exaFMM], [2.0])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_FILES([Makefile 2d/Makefile 2dp/Makefile 3d/Makefile 3dp/Makefile])
AC_CONFIG_HEADER(config.h)
AC_CONFIG_MACRO_DIR([m4])
#AC_CONFIG_SRCDIR([2d/fmm.cxx 2dp/fmm.cxx 3d/fmm.cxx 3dp/fmm.cxx])
AM_INIT_AUTOMAKE([foreign subdir-objects -Wall -Werror])
AM_MAINTAINER_MODE
AC_LANG([C++])
${CXXFLAGS=""}
${FCFLAGS=""}
AC_PROG_CXX([xlC_r xlC FCCpx FCC sxc++ CC icpc g++])
AC_PROG_FC([xlf90_r xlf90 frtpx frt sxf90 ftn ifort gfortran])
AC_PROG_RANLIB
AM_PROG_AR

# OpenMP
AC_ARG_ENABLE(openmp, [AC_HELP_STRING([--enable-openmp],[use OpenMP])], enable_openmp=$enableval, enable_openmp=no)
if test "$enable_openmp" = "yes"; then
   AC_DEFINE(EXAFMM_HAVE_OPENMP,1,[Use OpenMP])
   AX_OPENMP([], [AC_MSG_ERROR([don't know how to enable OpenMP for C++])])
fi
AM_CONDITIONAL(EXAFMM_HAVE_OPENMP, test "$enable_openmp" = "yes")

# Single precision
AC_ARG_ENABLE(single, [AC_HELP_STRING([--enable-single],[compile exaFMM in single precision])], ok=$enableval, ok=no)
AC_ARG_ENABLE(float, [AC_HELP_STRING([--enable-float],[synonym for --enable-single])], ok=$enableval)
if test "$ok" = "yes"; then
   AC_DEFINE(EXAFMM_SINGLE,1,[Define to compile in single precision.])
   PRECISION=s
else
   PRECISION=d
fi
AM_CONDITIONAL(EXAFMM_SINGLE, test "$ok" = "yes")
AC_SUBST(PRECISION)

# MPI
AC_ARG_ENABLE(mpi, [AC_HELP_STRING([--enable-mpi],[compile exaFMM MPI library])], enable_mpi=$enableval, enable_mpi=no)
if test "$enable_mpi" = "yes"; then
   AX_MPI([],[AC_MSG_ERROR([could not find mpi library for --enable-mpi])])
   AC_CHECK_PROG(MPIRUN, mpirun, mpirun)
   AC_SUBST(MPIRUN)
   CXX=$MPICXX
fi
AM_CONDITIONAL(EXAFMM_HAVE_MPI, test "$enable_mpi" = "yes")

# SSE3 (always defined)
AC_DEFINE(EXAFMM_HAVE_SSE3,1,[Define to enable SSE/SSE3 optimizations.])
if test "$PRECISION" != "d" -a "$PRECISION" != "s"; then
   AC_MSG_ERROR([SSE3 requires single or double precision])
fi
SSE3_CXXFLAGS=-msse3
AC_SUBST(SSE3_CXXFLAGS)

# AVX
AC_ARG_ENABLE(avx, [AC_HELP_STRING([--enable-avx],[enable AVX optimizations])], have_avx=$enableval, have_avx=no)
if test "$have_avx" = "yes"; then
   AC_DEFINE(EXAFMM_HAVE_AVX,1,[Define to enable AVX optimizations.])
   if test "$PRECISION" != "d" -a "$PRECISION" != "s"; then
      AC_MSG_ERROR([AVX requires single or double precision])
   fi
   AVX_CXXFLAGS=-mavx
fi
AM_CONDITIONAL(EXAFMM_HAVE_AVX, test "$have_avx" = "yes")
AC_SUBST(AVX_CXXFLAGS)

# AVX2
AC_ARG_ENABLE(avx2, [AC_HELP_STRING([--enable-avx2],[enable AVX2 optimizations])], have_avx2=$enableval, have_avx2=no)
if test "$have_avx2" = "yes"; then
   AC_DEFINE(EXAFMM_HAVE_AVX2,1,[Define to enable AVX2 optimizations.])
   if test "$PRECISION" != "d" -a "$PRECISION" != "s"; then
      AC_MSG_ERROR([AVX2 requires single or double precision])
   fi
   AVX2_CXXFLAGS=-mavx2
fi
AM_CONDITIONAL(EXAFMM_HAVE_AVX2, test "$have_avx2" = "yes")
AC_SUBST(AVX2_CXXFLAGS)

# MIC
AC_ARG_ENABLE(mic, [AC_HELP_STRING([--enable-mic],[enable MIC optimizations])], have_mic=$enableval, have_mic=no)
if test "$have_mic" = "yes"; then
   AC_DEFINE(EXAFMM_HAVE_MIC,1,[Define to enable MIC optimizations.])
   if test "$PRECISION" != "d" -a "$PRECISION" != "s"; then
      AC_MSG_ERROR([MIC requires single or double precision])
   fi
   MIC_CXXFLAGS=-mmic
fi
AM_CONDITIONAL(EXAFMM_HAVE_MIC, test "$have_mic" = "yes")
AC_SUBST(MIC_CXXFLAGS)

# AVX512
AC_ARG_ENABLE(avx512, [AC_HELP_STRING([--enable-avx512],[enable AVX512 optimizations])], have_avx512=$enableval, have_avx512=no)
if test "$have_avx512" = "yes"; then
   AC_DEFINE(EXAFMM_HAVE_AVX512,1,[Define to enable AVX512 optimizations.])
   if test "$PRECISION" != "d" -a "$PRECISION" != "s"; then
      AC_MSG_ERROR([AVX512 requires single or double precision])
   fi
   AVX512_CXXFLAGS=-xMIC-AVX512
fi
AM_CONDITIONAL(EXAFMM_HAVE_AVX512, test "$have_avx512" = "yes")
AC_SUBST(AVX512_CXXFLAGS)

# Neon
AC_ARG_ENABLE(neon, [AC_HELP_STRING([--enable-neon],[enable ARM NEON optimizations])], have_neon=$enableval, have_n\
eon=no)
if test "$have_neon" = "yes"; then
   AC_DEFINE(EXAFMM_HAVE_NEON,1,[Define to enable ARM NEON optimizations.])
   if test "$PRECISION" != "s"; then
      AC_MSG_ERROR([NEON requires single precision])
   fi
   NEON_CXXFLAGS=-mneon
fi
AM_CONDITIONAL(EXAFMM_HAVE_NEON, test "$have_neon" = "yes")
AC_SUBST(NEON_CXXFLAGS)

# Disable SIMD
AC_ARG_ENABLE(simd, [AC_HELP_STRING([--disable-simd],[disable SIMD optimizations])], use_simd=$enableval, use_simd=yes)
if test "$use_simd" = "no"; then
   AC_DEFINE(EXAFMM_USE_SIMD,0,[Disable SIMD optimizations.])
fi
AM_CONDITIONAL(EXAFMM_USE_SIMD, test "$use_simd" = "yes")

# Weights for partitioning
AC_ARG_ENABLE(weight, [AC_HELP_STRING([--enable-weight],[enable weight for partitioning])], use_weight=$enableval, use_weight=no)
if test "$use_weight" = "yes"; then
   AC_DEFINE(EXAFMM_USE_WEIGHT,1,[Enable weight for partitioning.])
fi
AM_CONDITIONAL(EXAFMM_USE_WEIGHT, test "$use_weight" = "yes")

# Kahan summation
AC_ARG_ENABLE(kahan, [AC_HELP_STRING([--enable-kahan],[enable Kahan summation])], use_kahan=$enableval, use_kahan=no)
if test "$use_kahan" = "yes"; then
   AC_DEFINE(EXAFMM_USE_KAHAN,1,[Enable Kahan summation.])
fi
AM_CONDITIONAL(EXAFMM_USE_KAHAN, test "$use_kahan" = "yes")

# PAPI
AC_ARG_ENABLE(papi, [AC_HELP_STRING([--enable-papi],[enable PAPI performance counter])], use_papi=$enableval, use_papi=no)
if test "$use_papi" = "yes"; then
   AC_DEFINE(EXAFMM_USE_PAPI,1,[Enable PAPI performance counter.])
fi
AM_CONDITIONAL(EXAFMM_USE_PAPI, test "$use_papi" = "yes")

# Trace
AC_ARG_ENABLE(trace, [AC_HELP_STRING([--enable-trace],[enable thread tracing])], use_trace=$enableval, use_trace=no)
if test "$use_trace" = "yes"; then
   AC_DEFINE(EXAFMM_USE_TRACE,1,[Enable thread tracing.])
fi
AM_CONDITIONAL(EXAFMM_USE_TRACE, test "$use_trace" = "yes")

# Count number of M2L and P2P kernel calls
AC_ARG_ENABLE(count-kernel, [AC_HELP_STRING([--enable-count-kernel],[count number of M2L and P2P kernel calls])], use_count_kernel=$enableval, use_count_kernel=no)
if test "$use_count_kernel" = "yes"; then
   AC_DEFINE(EXAFMM_COUNT_KERNEL,1,[Count number of M2L and P2P kernel calls.])
fi
AM_CONDITIONAL(EXAFMM_COUNT_KERNEL, test "$use_count_kernel" = "yes")

# Count interaction list
AC_ARG_ENABLE(count-list, [AC_HELP_STRING([--enable-count-list],[count interaction list per cell])], use_count_list=$enableval, use_count_list=no)
if test "$use_count_list" = "yes"; then
   AC_DEFINE(EXAFMM_COUNT_LIST,1,[Count interaction list per cell.])
fi
AM_CONDITIONAL(EXAFMM_COUNT_LIST, test "$use_count_list" = "yes")

# Assertion
AC_ARG_ENABLE(assert, [AC_HELP_STRING([--enable-assert],[enable assertion])], use_assert=$enableval, use_assert=no)
if test "$use_assert" = "yes"; then
   AC_DEFINE(EXAFMM_ASSERT,1,[Enable assertion.])
fi
AM_CONDITIONAL(EXAFMM_ASSERT, test "$use_assert" = "yes")

# Debug
AC_ARG_ENABLE(debug, [AC_HELP_STRING([--enable-debug],[compile with extra runtime checks for debugging])], enable_debug=$enableval, enable_debug=no)
if test "$enable_debug" = "yes"; then
   AC_DEFINE(EXAFMM_DEBUG,1,[Define to enable extra debugging options.])
fi
AM_CONDITIONAL(EXAFMM_DEBUG, test "$enable_debug" = "yes")

AX_COMPILER_VENDOR()
AX_COMPILER_FLAGS()
AC_LANG_PUSH([Fortran])
#AX_COMPILER_FLAGS()
AC_LANG_POP([Fortran])

AC_MSG_NOTICE([Vendor   : $ax_cv_cxx_compiler_vendor])
AC_MSG_NOTICE([CXX      : $CXX])
AC_MSG_NOTICE([CXXFLAGS : $COMPILER_CXXFLAGS])
AC_MSG_NOTICE([FC       : $FC])
AC_MSG_NOTICE([FCFLAGS  : $COMPILER_FCFLAGS])

AM_CONDITIONAL(EXAFMM_HAVE_BG, test "$ax_cv_cxx_compiler_vendor" = "ibm")
AM_CONDITIONAL(EXAFMM_HAVE_FX, test "$ax_cv_cxx_compiler_vendor" = "fujitsu")
AM_CONDITIONAL(EXAFMM_HAVE_SX, test "$ax_cv_cxx_compiler_vendor" = "nec")
AM_CONDITIONAL(EXAFMM_HAVE_CRAY, test "$CXX" = "CC")
AM_CONDITIONAL(EXAFMM_HAVE_INTEL, test "$ax_cv_cxx_compiler_vendor" = "intel")
AM_CONDITIONAL(EXAFMM_HAVE_CLANG, test "$ax_cv_cxx_compiler_vendor" = "clang")
AM_CONDITIONAL(EXAFMM_HAVE_GNU, test "$ax_cv_cxx_compiler_vendor" = "gnu")

AC_OUTPUT
