# ===========================================================================
#          http://www.gnu.org/software/autoconf-archive/ax_mpi.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_MPI([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro tries to find out how to compile programs that use MPI
#   (Message Passing Interface), a standard API for parallel process
#   communication (see http://www-unix.mcs.anl.gov/mpi/)
#
#   On success, it sets the MPICXX and MPIFC output variable
#   to the name of the MPI compiler, depending upon the current language.
#   (This may just be $CXX/$FC, but is more often something like
#   mpiCC/mpif90.) It also sets MPILIBS to any libraries that
#   are needed for linking MPI (e.g. -lmpi or -lfmpi, if a special
#   MPICXX/MPIFC was not found).
#
#   If you want to compile everything with MPI, you should use something
#   like this for C++:
#
#     if test -z "$CXX" && test -n "$MPICXX"; then
#       CXX="$MPICXX"
#     fi
#     AC_PROG_CXX
#     AX_MPI
#     CXX="$MPICXX"
#     LIBS="$MPILIBS $LIBS"
#
#   and similar for C++ (change all instances of CC to CXX), Fortran 77
#   (with F77 instead of CC) or Fortran (with FC instead of CC).
#
#   NOTE: The above assumes that you will use $CC (or whatever) for linking
#   as well as for compiling. (This is the default for automake and most
#   Makefiles.)
#
#   The user can force a particular library/compiler by setting the
#   MPICC/MPICXX/MPIF77/MPIFC and/or MPILIBS environment variables.
#
#   ACTION-IF-FOUND is a list of shell commands to run if an MPI library is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run if it is not
#   found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_MPI.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2008 Julian C. Cummings <cummings@cacr.caltech.edu>
#   Copyright (c) 2015 Rio Yokota <rioyokota@gmail.com>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

AC_DEFUN([AX_MPI], [
AC_PREREQ(2.50) dnl for AC_LANG_CASE
AC_REQUIRE([AC_PROG_CXX])
	AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
	AC_CHECK_PROGS(MPICXX, mpxlC_r mpxlC mpiFCCpx mpiFCC sxmpic++ mpiicpc mpicxx CC, $CXX)
	ax_mpi_save_CXX="$CXX"
	CXX="$MPICXX"
	AC_SUBST(MPICXX)
AC_REQUIRE([AC_PROG_FC])
	AC_ARG_VAR(MPIFC,[MPI Fortran compiler command])
	AC_CHECK_PROGS(MPIFC, mpxlf90_r mpxlf90 mpifrtpx mpifrt sxmpif03 ftn mpiifort mpif90, $FC)
	ax_mpi_save_FC="$FC"
	FC="$MPIFC"
	AC_SUBST(MPIFC)

if test x = x"$MPILIBS"; then
   AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])
   AC_CHECK_LIB(mpi, MPI_Init, [MPILIBS="-lmpi"])
   AC_CHECK_LIB(mpich, MPI_Init, [MPILIBS="-lmpich"])
   AC_CHECK_LIB(mpicxx, MPI_Init, [MPILIBS="-lmpicxx"])
   AC_CHECK_LIB(fmpi, MPI_Init, [MPILIBS="-lfmpi"])
   AC_CHECK_LIB(mpichf90, MPI_Init, [MPILIBS="-lmpichf90"])
fi
AC_CHECK_LIB(mpi_cxx, MPI_Init, [OPENMPILIBS="-lmpi_cxx"])
AC_SUBST(MPILIBS)
AC_SUBST(OPENMPILIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x = x"$MPILIBS"; then
        $2
        :
else
        ifelse([$1],,[AC_DEFINE(EXAFMM_HAVE_MPI,1,[Define if you have the MPI library.])],[$1])
        :
fi
])dnl AX_MPI
