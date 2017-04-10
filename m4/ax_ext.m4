# ===========================================================================
#          http://www.gnu.org/software/autoconf-archive/ax_ext.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_EXT
#
# DESCRIPTION
#
#   Find supported SIMD extensions by requesting cpuid. When an SIMD
#   extension is found, the -m"simdextensionname" is added to SIMD_FLAGS
#   (only if compilator support it) (ie : if "sse2" is available "-msse2" is
#   added to SIMD_FLAGS)
#
#   This macro calls:
#
#     AC_SUBST(SIMD_FLAGS)
#
# LICENSE
#
#   Copyright (c) 2008 Christophe Tournayre <turn3r@users.sourceforge.net>
#   Copyright (c) 2017 Tingyu Wang <twang66@gwu.edu>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

AC_DEFUN([AX_EXT],
[
  AC_REQUIRE([AC_CANONICAL_HOST])

  case $host_cpu in
    aarch64*)
      AC_CACHE_CHECK([whether NEON is supported], [ax_cv_have_neon_ext], [ax_cv_have_neon_ext=yes])
      if test "$ax_cv_have_neon_ext" = yes; then
        AX_CHECK_COMPILE_FLAG(-march=armv8-a+simd, [SIMD_FLAGS="-march=armv8-a+simd"], [ax_cv_have_neon_ext=no])
      fi
      ;;

    arm*)
      AC_CACHE_CHECK([whether NEON is supported], [ax_cv_have_neon_ext], [ax_cv_have_neon_ext=yes])
      if test "$ax_cv_have_neon_ext" = yes; then
        AX_CHECK_COMPILE_FLAG(-mfpu=neon, [SIMD_FLAGS="-mfpu=neon"], [ax_cv_have_neon_ext=no])
      fi
      ;;

    i[[3456]]86*|x86_64*|amd64*)
      AC_CACHE_CHECK([whether sse3 is supported], [ax_cv_have_sse3_ext], [ax_cv_have_sse3_ext=yes])
      if test "$ax_cv_have_sse3_ext" = yes; then
        AX_CHECK_COMPILE_FLAG(-msse3, [SIMD_FLAGS="-msse3"], [ax_cv_have_sse3_ext=no])
      fi

      AC_CACHE_CHECK([whether avx is supported], [ax_cv_have_avx_ext], [ax_cv_have_avx_ext=yes])
      if test "$ax_cv_have_avx_ext" = yes; then
        AX_CHECK_COMPILE_FLAG(-mavx, [SIMD_FLAGS="-mavx"], [ax_cv_have_avx_ext=no])
      fi

      AC_CACHE_CHECK([whether avx2 is supported], [ax_cv_have_avx2_ext], [ax_cv_have_avx2_ext=yes])
      if test "$ax_cv_have_avx2_ext" = yes; then
        AX_CHECK_COMPILE_FLAG(-mavx2, [SIMD_FLAGS="-mavx2"], [ax_cv_have_avx2_ext=no])
      fi
      ;;
  esac

  AC_SUBST(SIMD_FLAGS)
])
