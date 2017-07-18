AC_PREREQ([2.69])
AC_INIT([exaFMM], [2.0])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_FILES([Makefile 2d/Makefile 3d/Makefile 3dp/Makefile 3dp_gpu/Makefile 3dp_gpu_mpi/Makefile 3dp_simd/Makefile 3dp_simd_mpi/Makefile])
AC_CONFIG_HEADER(config.h)
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign subdir-objects -Wall -Werror])
AM_MAINTAINER_MODE
AC_LANG([C++])
#: ${CXXFLAGS=""}
#: ${FCFLAGS=""}
AC_PROG_CXX([CC icpc g++])
AC_PROG_FC([ftn ifort gfortran])
AC_PROG_RANLIB
AM_PROG_AR

# OpenMP
AC_ARG_ENABLE([openmp],
              [AS_HELP_STRING([--disable-openmp],[disable OpenMP])],
              [],
              [enable_openmp=yes])
if test "$enable_openmp" = "yes"; then
   AC_DEFINE(EXAFMM_HAVE_OPENMP,1,[Use OpenMP])
   AX_OPENMP([], [AC_MSG_ERROR([don't know how to enable OpenMP for C++])])
fi
AM_CONDITIONAL(EXAFMM_HAVE_OPENMP, test "$enable_openmp" = "yes")

# MPI
AC_ARG_ENABLE([mpi],
              [AS_HELP_STRING([--disable-mpi],[do not compile with MPI])],
              [],
              [enable_mpi=yes])
if test "$enable_mpi" = "yes"; then
   AX_MPI([],[AC_MSG_ERROR([could not find mpi library for --enable-mpi])])
   AC_CHECK_PROG(MPIRUN, mpirun, mpirun)
   AC_SUBST(MPIRUN)
   CXX=$MPICXX
fi
AM_CONDITIONAL(EXAFMM_HAVE_MPI, test "$enable_mpi" = "yes")

# SIMD extensions
AC_ARG_ENABLE([simd],
              [AS_HELP_STRING([--disable-simd],[disable SIMD])],
              [enable_simd=$enableval],
              [enable_simd=yes])
if test "$enable_simd" = "yes"; then
   AX_EXT
fi

# CUDA
AC_ARG_ENABLE([cuda],
              [AS_HELP_STRING([--disable-cuda],[do not compile with CUDA])],
              [],
              [enable_cuda=yes])
if test "$enable_cuda" = "yes"; then
   AC_DEFINE(EXAFMM_HAVE_CUDA,1,[Use CUDA])
   AX_CHECK_CUDA([],[AC_MSG_ERROR([could not find CUDA runtime])])
fi
AM_CONDITIONAL(EXAFMM_HAVE_CUDA, test "$enable_cuda" = "yes")

# Single precision
AC_ARG_ENABLE([single],
              [AS_HELP_STRING([--enable-single],[compile in single precision])],
              [enable_single=$enableval],
              [enable_single=no])
AC_ARG_ENABLE([float],
              [AS_HELP_STRING([--enable-float],[synonym for --enable-single])],
              [enable_single=$enableval])
if test "$enable_single" = "yes"; then
   AC_DEFINE(EXAFMM_SINGLE,1,[Define to compile in single precision.])
   PRECISION=s
else
   PRECISION=d
fi
AM_CONDITIONAL(EXAFMM_SINGLE, test "$enable_single" = "yes")
AC_SUBST(PRECISION)

# Assertion
AC_ARG_ENABLE([assert],
              [AS_HELP_STRING([--disable-assert],[disable assertion])],
              [],
              [enable_assert=yes])
if test "$enable_assert" = "yes"; then
   AC_DEFINE(EXAFMM_ASSERT,1,[Enable assertion.])
fi
AM_CONDITIONAL(EXAFMM_ASSERT, test "$use_assert" = "yes")

# Debug
AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],[enable debugging])],
              [enable_debug=$enableval],
              [enable_debug=no])
if test "$enable_debug" = "yes"; then
   AC_DEFINE(EXAFMM_DEBUG,1,[Enable debugging options.])
fi
AM_CONDITIONAL(EXAFMM_DEBUG, test "$enable_debug" = "yes")

AX_COMPILER_VENDOR()
AX_COMPILER_FLAGS()

AC_MSG_NOTICE([OpenMP   : $enable_openmp])
AC_MSG_NOTICE([MPI      : $enable_mpi])
AC_MSG_NOTICE([CUDA     : $enable_cuda])
AC_MSG_NOTICE([FP32     : $enable_single])
AC_MSG_NOTICE([SIMD     : $SIMD_FLAGS])
AC_MSG_NOTICE([Assert   : $enable_assert])
AC_MSG_NOTICE([Debug    : $enable_debug])
AC_MSG_NOTICE([Vendor   : $ax_cv_cxx_compiler_vendor])
AC_MSG_NOTICE([CXX      : $CXX])
AC_MSG_NOTICE([FC       : $FC])
AC_MSG_NOTICE([CXXFLAGS : $CXXFLAGS])
AC_MSG_NOTICE([FCFLAGS  : $FCFLAGS])

AM_CONDITIONAL(EXAFMM_HAVE_CRAY, test "$CXX" = "CC")
AM_CONDITIONAL(EXAFMM_HAVE_INTEL, test "$ax_cv_cxx_compiler_vendor" = "intel")
AM_CONDITIONAL(EXAFMM_HAVE_CLANG, test "$ax_cv_cxx_compiler_vendor" = "clang")
AM_CONDITIONAL(EXAFMM_HAVE_GNU, test "$ax_cv_cxx_compiler_vendor" = "gnu")

AC_OUTPUT
