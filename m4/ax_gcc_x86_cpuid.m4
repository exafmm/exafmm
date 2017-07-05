#serial 10

AC_DEFUN([AX_GCC_X86_CPUID],
[AX_GCC_X86_CPUID_COUNT($1, 0)
])

AC_DEFUN([AX_GCC_X86_CPUID_COUNT],
[AC_REQUIRE([AC_PROG_CC])
AC_LANG_PUSH([C])
AC_CACHE_CHECK(for x86 cpuid $1 output, ax_cv_gcc_x86_cpuid_$1,
 [AC_RUN_IFELSE([AC_LANG_PROGRAM([#include <stdio.h>], [
     int op = $1, level = $2, eax, ebx, ecx, edx;
     FILE *f;
      __asm__ __volatile__ ("xchg %%ebx, %1\n"
        "cpuid\n"
        "xchg %%ebx, %1\n"
        : "=a" (eax), "=r" (ebx), "=c" (ecx), "=d" (edx)
        : "a" (op), "2" (level));

     f = fopen("conftest_cpuid", "w"); if (!f) return 1;
     fprintf(f, "%x:%x:%x:%x\n", eax, ebx, ecx, edx);
     fclose(f);
     return 0;
])],
     [ax_cv_gcc_x86_cpuid_$1=`cat conftest_cpuid`; rm -f conftest_cpuid],
     [ax_cv_gcc_x86_cpuid_$1=unknown; rm -f conftest_cpuid],
     [ax_cv_gcc_x86_cpuid_$1=unknown])])
AC_LANG_POP([C])
])
