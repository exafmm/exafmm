#ifndef namespace_h
#define namespace_h
#ifndef EXAFMM_NAMESPACE
#if EXAFMM_LAPLACE
#define EXAFMM_NAMESPACE exafmm_laplace
#elif EXAFMM_HELMHOLTZ
#define EXAFMM_NAMESPACE exafmm_helmholtz
#elif EXAFMM_BIOTSAVART
#define EXAFMM_NAMESPACE exafmm_biotsavart
#else
#error Please define macro EXAFMM_LAPLACE or EXAFMM_HELMHOLTZ or EXAFMM_BIOTSAVART
#endif
#endif
#endif
