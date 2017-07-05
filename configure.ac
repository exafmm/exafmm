AC_PREREQ([2.69])
AC_INIT([exaFMM], [2.0])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_FILES([Makefile 2d/Makefile 3d/Makefile 3dp/Makefile 3dp_gpu/Makefile 3dp_gpu_mpi/Makefile 3dp_simd/Makefile 3dp_simd_mpi/Makefile])
AC_CONFIG_HEADER(config.h)
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign subdir-objects -Wall -Werror])
AM_MAINTAINER_MODE
AC_LANG([C++])
${CXXFLAGS=""}
${FCFLAGS=""}
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

# Single precision
AC_ARG_ENABLE([single],
              [AS_HELP_STRING([--enable-single],[compile in single precision])],
              [ok=$enableval],
              [ok=no])
AC_ARG_ENABLE([float],
              [AS_HELP_STRING([--enable-float],[synonym for --enable-single])],
              [ok=$enableval])
if test "$ok" = "yes"; then
   AC_DEFINE(EXAFMM_SINGLE,1,[Define to compile in single precision.])
   PRECISION=s
else
   PRECISION=d
fi
AM_CONDITIONAL(EXAFMM_SINGLE, test "$ok" = "yes")
AC_SUBST(PRECISION)

# SIMD extensions supported by compiler
AC_ARG_ENABLE([simd-ext],
              [AS_HELP_STRING([--disable-simd-ext],[disable SIMD])],
              [use_simd_ext=$enableval],
              [use_simd_ext=yes])
if test "$use_simd_ext" = "yes"; then
   AX_EXT
fi

# Assertion
AC_ARG_ENABLE([assert],
              [AS_HELP_STRING([--enable-assert],[enable assertion])],
              [use_assert=$enableval],
              [use_assert=no])
if test "$use_assert" = "yes"; then
   AC_DEFINE(EXAFMM_ASSERT,1,[Enable assertion.])
fi
AM_CONDITIONAL(EXAFMM_ASSERT, test "$use_assert" = "yes")

# Debug
AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],[compile in debugging mode])],
              [enable_debug=$enableval],
              [enable_debug=no])
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
AC_MSG_NOTICE([FC       : $FC])

AM_CONDITIONAL(EXAFMM_HAVE_CRAY, test "$CXX" = "CC")
AM_CONDITIONAL(EXAFMM_HAVE_INTEL, test "$ax_cv_cxx_compiler_vendor" = "intel")
AM_CONDITIONAL(EXAFMM_HAVE_CLANG, test "$ax_cv_cxx_compiler_vendor" = "clang")
AM_CONDITIONAL(EXAFMM_HAVE_GNU, test "$ax_cv_cxx_compiler_vendor" = "gnu")

AC_OUTPUT
